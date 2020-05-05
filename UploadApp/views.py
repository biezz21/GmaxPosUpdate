from django.shortcuts import render
from .models import V1Model, V2Model
from django.core.files.storage import FileSystemStorage
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
from django.conf import settings
import subprocess

# Create your views here.

def UploadView(request):
    if request.method == 'POST' and request.FILES['vcf_file']:
        vcf_file = request.FILES['vcf_file']
        fs = FileSystemStorage()
        filename = fs.save(vcf_file.name, vcf_file)
        uploaded_file_url = fs.url(filename)
        refversion = request.POST.get('Refversion')

        # V4
        # ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm4.4PTR.genome_main.fna'))
        # f           = ref_file.read()
        # ref         = f.split('>')

        # dicfa_v4    = {}
        # for x in ref[1:]:
        #     title           = x.split('\n')[0].split()[0].split('.')[-1]
        #     seq             = ''.join(x.split('\n')[1:])
        #     dicfa_v4[title] = seq

        # Select Ref Version
        if refversion == 'Glycine max accession Williams genome assembly v1.0':
            ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm1.FCtY.genome_main.fna'))
            df = pd.DataFrame(list(V1Model.objects.all().values()), columns=['Pos_v4','Pos_pre','Pos_pre_Info'])
            ix = list(df['Pos_pre'])
        elif refversion == 'Glycine max accession Williams 82 genome assembly v2.0':
            ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm2.DTC4.genome_main.fna'))
            df = pd.DataFrame(list(V1Model.objects.all().values()), columns=['Pos_v4','Pos_pre','Pos_pre_Info'])
            ix = list(df['Pos_pre'])
            
        f           = ref_file.read()
        ref         = f.split('>')

        dicfa_pre   = {}
        for x in ref[1:]:
            title            = x.split('\n')[0].split()[0].split('.')[-1]
            seq              = ''.join(x.split('\n')[1:])
            dicfa_pre[title] = seq

        vcf_mat = pd.read_csv('%s'%os.path.join(settings.MEDIA_ROOT, filename), compression='gzip',
                comment='#', sep='\t', header=None)

        vcf_mat = vcf_mat.dropna(axis=0)

        def NormalChr(x):
            if x.split('_')[0] == 'SCAFFOLD':
                return x.split('_')[0].lower()+ '_' + x.split('_')[1]
            else:
                try:
                    return 'Gm' + x.split('Chr')[1]
                except IndexError:
                    return x

        vcf_mat[0] = [NormalChr(x.split('.')[-1]) for x in vcf_mat[0]]

        # Filtering INDEL

        a = np.array([len(x.split(',')[0]) for x in vcf_mat[4]])
        b = np.array([len(x) for x in vcf_mat[3]])

        m_alt = ( a < 2 )
        m_ref = (b < 2)
        m_indel = ['INDEL' not in x.split(';') for x in vcf_mat[7]]

        m = m_alt & m_ref & m_indel
        vcf_mat_snp = vcf_mat[m]

        vcf_mat_snp['ix'] = vcf_mat_snp.apply(lambda x : x[0]+'-'+str(x[1]), axis=1)
        vcf_mat_snp_ix = vcf_mat_snp.set_index('ix')

        common = set(ix).intersection(list(vcf_mat_snp['ix']))
        Noncommon = set(list(vcf_mat_snp['ix'])).difference(ix)

        vcf_mat_common = vcf_mat_snp_ix.loc[common]
        df_set = df.set_index('Pos_pre').loc[common]

        chrom = [x.split('-')[0] for x in df_set['Pos_v4']]
        pos = [x.split('-')[1] for x in df_set['Pos_v4']]

        vcf_mat_common[0] = chrom
        vcf_mat_common[1] = pos

        def ReversePosUpdated(x):
            dic_cg = {'A':'T','G':'C','T':'A', 'C':'G'}
            
            
            if df.loc[x.index]['Pos_pre_Info'] == False:
                
                FixREF.append(dic_cg[x[3]])
                FixALT.append(','.join([dic_cg[i] for i in x[4].split(',')]))
                
                
                
            else:
                FixREF.append(x[3])
                FixALT.append(x[4])

        global FixREF, FixALT
        FixREF, FixALT = [], []

        abc = vcf_mat_common.apply(lambda x : ReversePosUpdated(x), axis=1)

        vcf_mat_common[3] = FixREF
        vcf_mat_common[4] = FixALT

        
        
        
        
        obj_output = {
            'uploaded_file_url': uploaded_file_url,
            'refversion': refversion,
        }

        fs.delete(vcf_file.name)
        return render(request, 'download.html', obj_output)

    return render(request, 'upload.html')