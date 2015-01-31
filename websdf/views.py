from django.shortcuts import render, render_to_response
from .forms import UploadFileForm
from .calculations import read_sdf

import pandas as pd

def home(request):
    '''
    Define home page view. It is loaded by urls.py.
    '''
    
    return render(request, 'index.html', {'form':UploadFileForm})

def upload_file(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            df = read_sdf(request.FILES['file'])
            cols = list(df.columns)
            
            rows = []
            for ix, row in df.iterrows():
            #for cell in row:
            #    print cell
                rows.append(list(row))
            print rows
            return render(request, 'table.html', {'rows':rows, 'cols':cols})
        

