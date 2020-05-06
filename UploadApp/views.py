from django.shortcuts import render
from .models import V1Model, V2Model
from django.core.files.storage import FileSystemStorage
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os, sys
from django.conf import settings
import gzip
import shutil

# Create your views here.

def UploadView(request):
    if request.method == 'POST' and request.FILES['vcf_file']:
        vcf_file = request.FILES['vcf_file']
        fs = FileSystemStorage()
        filename = fs.save(vcf_file.name, vcf_file)
        refversion = request.POST.get('Refversion')

        # Select Ref Version
        if refversion == 'Glycine max accession Williams genome assembly v1.0':
            ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm1.FCtY.genome_main.fna'))
            df = pd.DataFrame(list(V1Model.objects.all().values()))
            ix = list(df['Pos_v1'])
        elif refversion == 'Glycine max accession Williams 82 genome assembly v2.0':
            ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm2.DTC4.genome_main.fna'))
            df = pd.DataFrame(list(V2Model.objects.all().values()))
            ix = list(df['Pos_v2'])
            
        f           = ref_file.read()
        ref         = f.split('>')

        dicfa_pre   = {}
        for x in ref[1:]:
            title            = x.split('\n')[0].split()[0].split('.')[-1]
            seq              = ''.join(x.split('\n')[1:])
            dicfa_pre[title] = seq

        df.columns = ['id', 'Pos_v4', 'Pos_pre', 'Pos_pre_Info']

        vcf_mat = pd.read_csv('%s'%os.path.join(settings.MEDIA_ROOT, filename), compression='gzip',
                comment='#', sep='\t', header=None)

        import subprocess

        c = subprocess.check_output(["sh","%s"%os.path.join(settings.MEDIA_ROOT, 'columns.sh'),"%s"%os.path.join(settings.MEDIA_ROOT, filename)], universal_newlines=True)
        vcf_mat.columns = c.split('\n')[:-1][0].split('\t')

        vcf_mat = vcf_mat.dropna(axis=0)

        # vcf_mat.to_csv('%s'% os.path.join(settings.MEDIA_ROOT, 'test.csv'))

        # os.system('less %s | grep "#" > %s'%(os.path.join(settings.MEDIA_ROOT, filename), os.path.join(settings.MEDIA_ROOT, 'test.txt')))

        def NormalChr(x):
            if x.split('_')[0] == 'SCAFFOLD':
                return x.split('_')[0].lower()+ '_' + x.split('_')[1]
            else:
                try:
                    return 'Gm' + x.split('Chr')[1]
                except IndexError:
                    return x

        vcf_mat['#CHROM'] = [NormalChr(x.split('.')[-1]) for x in vcf_mat['#CHROM']]

        # # Filtering INDEL

        a = np.array([len(x.split(',')[0]) for x in vcf_mat['ALT']])
        b = np.array([len(x) for x in vcf_mat['REF']])

        m_alt = ( a < 2 )
        m_ref = (b < 2)
        m_indel = ['INDEL' not in x.split(';') for x in vcf_mat['INFO']]

        m = m_alt & m_ref & m_indel
        vcf_mat_snp = vcf_mat[m]

        vcf_mat_snp['ix'] = vcf_mat_snp.apply(lambda x : x['#CHROM']+'-'+str(x['POS']), axis=1)
        vcf_mat_snp_ix = vcf_mat_snp.set_index('ix')
        vcf_mat_snp_ix_info = vcf_mat_snp_ix[vcf_mat_snp_ix.columns[0:9]]
        vcf_mat_snp_ix_samples = vcf_mat_snp_ix[vcf_mat_snp_ix.columns[9:]]
        

        common = set(ix).intersection(list(vcf_mat_snp['ix']))
        Noncommon = set(list(vcf_mat_snp['ix'])).difference(ix)

        vcf_mat_common = vcf_mat_snp_ix_info.loc[common]
        df_set = df.set_index('Pos_pre').loc[common]

        chrom = [x.split('-')[0] for x in df_set['Pos_v4']]
        pos = [x.split('-')[1] for x in df_set['Pos_v4']]

        vcf_mat_common['#CHROM'] = chrom
        vcf_mat_common['POS'] = pos

        def ReversePosUpdated(x):
            dic_cg = {'A':'T','G':'C','T':'A', 'C':'G'}
            
            
            if df_set.loc[x.name]['Pos_pre_Info'] == False:
                
                FixREF.append(dic_cg[x['REF']])
                FixALT.append(','.join([dic_cg[i] for i in x['ALT'].split(',')]))
                
                
                
            else:
                FixREF.append(x['REF'])
                FixALT.append(x['ALT'])

        global FixREF, FixALT
        FixREF, FixALT = [], []

        progress = vcf_mat_common.apply(lambda x : ReversePosUpdated(x), axis=1)

        vcf_mat_common['REF'] = FixREF
        vcf_mat_common['ATL'] = FixALT

        vcf_common = pd.merge(vcf_mat_common, vcf_mat_snp_ix_samples.loc[common], left_index=True, right_index=True, how='left')

        # vcf_common.to_csv('%s'%os.path.join(settings.MEDIA_ROOT, 'test.csv'), index=False)

        vcf_mat_None = vcf_mat_snp_ix_info.loc[Noncommon]

        # vcf_mat_None.to_csv('%s'%os.path.join(settings.MEDIA_ROOT, 'test.csv'), index=False)

        def MakeSeq(matrix):
    
            try:
                pos = matrix['POS']-1
                
                if len(dicfa_pre[matrix['#CHROM']][pos-500:pos+500]) == 1000:
                    
                    Target_seq     = dicfa_pre[matrix['#CHROM']][pos-500:pos+500]
                    Seq_file       = SeqRecord(Seq(Target_seq), id="%s-%s"%(matrix['#CHROM'],str(matrix['POS'])), 
                                            name='', description = '')
                    
                    seq_all.append(Seq_file)
                    
                    return "Done"
                else:
                    if pos < 500 and len(dicfa_pre[matrix['#CHROM']][pos:pos+1000]) == 1000:
                        Target_seq     = dicfa_pre[matrix['#CHROM']][pos:pos+1000]
                        Seq_file       = SeqRecord(Seq(Target_seq), id="%s-%s"%(matrix['#CHROM'],str(matrix['POS'])), 
                                            name='', description = '')
                        
                        seq_all.append(Seq_file)
                        
                        return "x<500"
                    
                    elif len(dicfa_pre[matrix['#CHROM']][pos-999:pos+1]) == 1000:
                        Target_seq   = dicfa_pre[matrix['#CHROM']][pos-999:pos+1]
                        Seq_file       = SeqRecord(Seq(Target_seq), id="%s-%s"%(matrix['#CHROM'],str(matrix['POS'])), 
                                            name='', description = '')
                        
                        seq_all.append(Seq_file)
                        
                        return "500<x"
                    else:
                        return "%s-%s"%(matrix['#CHROM'],str(matrix['POS']))
                    
                    
            except KeyError:
                Err.append(["%s-%s"%(matrix['#CHROM'],str(matrix['POS']))])
                
                
                
                
        global Err, seq_all
        Err, seq_all = [], []

        PositionInfo = vcf_mat_None.apply(lambda x: MakeSeq(x), axis=1)
        SeqIO.write(seq_all, '%s'%os.path.join(settings.MEDIA_ROOT, 'sample.fasta'), 'fasta')

        os.system('bwa mem -t 30 -A 100 -k 1000 %s %s > %s'%(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm4.4PTR.genome_main.fna'), os.path.join(settings.MEDIA_ROOT, 'sample.fasta'), os.path.join(settings.MEDIA_ROOT, 'sample.sam')))
        df_bwa = pd.read_csv('%s'%os.path.join(settings.MEDIA_ROOT, 'sample.sam'), sep='\t', header=None, comment='@', names=range(16))
        df_bwa = df_bwa[df_bwa.columns[[0,2,3,5]]]
        df_bwa.columns = ['id', '#CHROM', 'POS', 'align_length']
        m = df_bwa['align_length'] == '1000M'
        df_bwa_matched = df_bwa[m]

        if df_bwa_matched.shape[0] == 0:
            pass
        else:
            # V4
            ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm4.4PTR.genome_main.fna'))
            f           = ref_file.read()
            ref         = f.split('>')

            dicfa_v4    = {}
            for x in ref[1:]:
                title           = x.split('\n')[0].split()[0].split('.')[-1]
                seq             = ''.join(x.split('\n')[1:])
                dicfa_v4[title] = seq

            vcf_mat_None['PositionInfo'] = PositionInfo
            bwa_id = list(df_bwa_matched['id'])
            vcf_updated = vcf_mat_None.loc[bwa_id]

            bwa_CHROM = [x.split('.')[-1] for x in df_bwa_matched['#CHROM']]
            df_bwa_matched['#CHROM'] = bwa_CHROM

            def PositionUpdate(x):

                if vcf_updated.loc[x['id']]['PositionInfo'] == 'x<500':

                    if dicfa_v4[x['#CHROM']][x['POS']-1] == vcf_updated.loc[x['id']]['REF']:

                        ForwardInfo.append(True)
                        return x['POS']

                    else:

                        ForwardInfo.append(False)
                        return x['POS'] + 999

                elif vcf_updated.loc[x['id']]['PositionInfo'] == '500<x':

                    if dicfa_v4[x['#CHROM']][x['POS']+999-1] == vcf_updated.loc[x['id']]['REF']:

                        ForwardInfo.append(True)
                        return x['POS'] + 999

                    else:

                        ForwardInfo.append(False)
                        return x['POS']

                else:

                    if dicfa_v4[x['#CHROM']][x['POS']+500-1] == vcf_updated.loc[x['id']]['REF']:

                        ForwardInfo.append(True)
                        return x['POS'] + 500

                    else:

                        ForwardInfo.append(False)
                        return x['POS'] + 500 -1

            global ForwardInfo
            ForwardInfo = []
            POS = df_bwa_matched.apply(lambda x : PositionUpdate(x), axis=1)
            CHROM = list(df_bwa_matched['#CHROM'])

            vcf_updated = vcf_updated[vcf_updated.columns[:-1]]

            vcf_updated['#CHROM'] = CHROM
            vcf_updated['POS'] = POS
            vcf_updated['strand'] = ForwardInfo
            Pos_v4 = list(vcf_updated.index)


            def ReversePosUpdated(x):
                dic_cg = {'A':'T','G':'C','T':'A', 'C':'G'}
                
                
                if x['strand'] == False:
                    
                    FixREF2.append(dic_cg[x['REF']])
                    FixALT2.append(','.join([dic_cg[i] for i in x['ALT'].split(',')]))
                    
                    
                    
                else:
                    FixREF2.append(x['REF'])
                    FixALT2.append(x['ALT'])

            global FixREF2, FixALT2
            FixREF2, FixALT2 = [], []

            progress = vcf_updated.apply(lambda x : ReversePosUpdated(x), axis=1)

            vcf_updated['REF'] = FixREF2
            vcf_updated['ALT'] = FixALT2
            vcf_updated = vcf_updated[vcf_updated.columns[:-1]]
            vcf_non = pd.merge(vcf_updated, vcf_mat_snp_ix_samples.loc[bwa_id], left_index=True, right_index=True, how='left')
            # vcf_non.to_csv('%s'%os.path.join(settings.MEDIA_ROOT, 'test.csv'), index=False)
            chr_pos = [x+'-'+str(pos) for x,y in zip(CHROM, pos)]
            if refversion == 'Glycine max accession Williams genome assembly v1.0':
                for x,y,z in zip(Pos_v4, chr_pos, ForwardInfo):
                    model = V1Model(
                        Pos_v4 = x,
                        Pos_v1 = y,
                        Pos_v1_Info = z,
                    )
                    model.save()
            elif refversion == 'Glycine max accession Williams 82 genome assembly v2.0':
                for x,y,z in zip(Pos_v4, chr_pos, ForwardInfo):
                    model = V1Model(
                        Pos_v4 = x,
                        Pos_v2 = y,
                        Pos_v2_Info = z,
                    )
                    model.save()




        try:
            df_result = pd.concat([vcf_common, vcf_non])
        except UnboundLocalError:
            df_result = vcf_common

        # df_result.to_csv('%s'%os.path.join(settings.MEDIA_ROOT, 'test.csv'), index=False)

        header = '''##reference=glyma.Wm82.gnm4.4PTR.genome_main.fna\n'''
        output_VCF = '%s'%os.path.join(settings.MEDIA_ROOT, 'output.vcf')
        with open(output_VCF, 'w') as vcf:
            vcf.write(header)

        df_result.to_csv(output_VCF, sep="\t", mode='a', index=False)

        # with open('%s'%os.path.join(settings.MEDIA_ROOT, 'output.vcf'), 'rb') as f_in:
        #     with gzip.open('%s'%os.path.join(settings.MEDIA_ROOT, 'output.vcf.gz'), 'wb') as f_out:
        #         shutil.copyfileobj(f_in, f_out)

        
        uploaded_file_url = fs.url('output.vcf')
            


        
        obj_output = {
            'uploaded_file_url': uploaded_file_url,
            }

        fs.delete(vcf_file.name)
        return render(request, 'upload.html', obj_output)

    return render(request, 'upload.html')