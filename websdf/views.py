from django.shortcuts import render, render_to_response, redirect
from .forms import UploadFileForm
from .calculations import read_sdf

import pandas as pd

def home(request):
    '''
    Define home page view. It is loaded by urls.py.
    '''
    if 'error' in request.GET:
        return render(request, 'index.html', {'form':UploadFileForm, 'error':request.GET['error']})
    else:
        return render(request, 'index.html', {'form':UploadFileForm})

def upload_file(request):
    '''
    View that returns web page with table with SDF contents
    '''
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid() and str(request.FILES['file']).endswith('.sdf'):
            checks = request.POST.getlist('checks')
            df = read_sdf(request.FILES['file'], checks)
                
            # Get table columns
            cols = list(df.columns)
            
            # Extract DataFrame rows
            rows = []
            for ix, row in df.iterrows():
                rows.append(zip(cols,list(row)))
            print checks
            #print rows
            return render(request, 'table.html', {'rows':rows, 'cols':cols})
        else:
            error = 'Please select a valid SDF file'
        
    else:
        error = 'Please select a valid SDF file'
    return redirect('/?error=%s' %error)
        

