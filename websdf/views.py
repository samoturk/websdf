from django.shortcuts import render, redirect
from websdf.calculations import read_sdf, read_smi
from websdf.settings import PAGE_URL

def home(request):
    '''
    Define home page view. It is loaded by urls.py.
    '''
    if 'error' in request.GET:
        return render(request, 'index.html', {'PAGE_URL':PAGE_URL, 
                                              'error':request.GET['error']})
    else:
        return render(request, 'index.html', {'PAGE_URL':PAGE_URL})

def upload_file(request):
    '''
    View that returns web page with table with SDF contents
    '''
    if request.method == 'POST':
        if str(request.FILES['file']).endswith('.sdf'):
            checks = request.POST.getlist('checks')
            df = read_sdf(request.FILES['file'], checks)
            # Get table columns
            cols = list(df.columns)
           # Extract DataFrame rows
            rows = []
            for ix, row in df.iterrows():
                rows.append(zip(cols,list(row)))
            return render(request, 'table.html', {'PAGE_URL':PAGE_URL,
                                                  'rows':rows, 'cols':cols})
            
        elif (str(request.FILES['file']).endswith('.smi') or str(request.FILES[
            'file']).endswith('.ism')):
            checks = request.POST.getlist('checks')            
            df = read_smi(request.FILES['file'], checks)            
            # Get table columns
            cols = list(df.columns)
            # Extract DataFrame rows
            rows = []
            for ix, row in df.iterrows():
                rows.append(zip(cols,list(row)))
            return render(request, 'table.html', {'PAGE_URL':PAGE_URL,
                                                  'rows':rows, 'cols':cols})
        else:
            error = 'Please select a valid SDF/SMI file'
        
    else:
        error = 'Please select a valid SDF/SMI file'
    return redirect(PAGE_URL + '/?error=%s' %error)
        

