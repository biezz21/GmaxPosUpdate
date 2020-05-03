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
        elif refversion == 'Glycine max accession Williams 82 genome assembly v2.0':
            ref_file    = open(os.path.join(settings.MEDIA_ROOT, 'glyma.Wm82.gnm2.DTC4.genome_main.fna'))
            
        f           = ref_file.read()
        ref         = f.split('>')

        dicfa_pre   = {}
        for x in ref[1:]:
            title            = x.split('\n')[0].split()[0].split('.')[-1]
            seq              = ''.join(x.split('\n')[1:])
            dicfa_pre[title] = seq

        vcf_mat = pd.read_csv('%s'%os.path.join(settings.MEDIA_ROOT, filename), compression='gzip',
                comment='#', sep='\t', header=None)
            
        obj_output = {
            'uploaded_file_url': uploaded_file_url,
            'refversion': refversion,
        }

        fs.delete(vcf_file.name)
        return render(request, 'download.html', obj_output)

    return render(request, 'upload.html')