from django.shortcuts import render
from .models import V1Model, V2Model
from django.core.files.storage import FileSystemStorage

# Create your views here.

def UploadView(request):
    if request.method == 'POST' and request.FILES['vcf_file']:
        vcf_file = request.FILES['vcf_file']
        fs = FileSystemStorage()
        filename = fs.save(vcf_file.name, vcf_file)
        uploaded_file_url = fs.url(filename)
        refversion = request.POST.get('Refversion')

        obj_output = {
            'uploaded_file_url': uploaded_file_url,
            'refversion': refversion,
        }

        fs.delete(vcf_file.name)
        return render(request, 'download.html', obj_output)

    return render(request, 'upload.html')