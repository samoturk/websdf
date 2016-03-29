from django.shortcuts import render, redirect
from django.http import JsonResponse, HttpResponseNotAllowed, HttpResponseBadRequest
from websdf.calculations import read_sdf, read_smi, read_mol, read_mol2, read_smi_string
from websdf.settings import PAGE_URL
from rdkit.Chem import MolFromSmiles
try:
    # Import extra, proprietary functions
    from websdf.extra import extra
except:
    extra = None
    
def home(request):
    """
    Define home page view. It is loaded by urls.py.
    """
    if 'error' in request.GET:
        return render(request, 'index.html', {'PAGE_URL':PAGE_URL, 
                                              'error':request.GET['error']})
    else:
        return render(request, 'index.html', {'PAGE_URL':PAGE_URL})

def upload_file(request):
    """
    View that returns web page with table with SDF contents
    """
    if request.method == 'POST':
        filename = str(request.FILES['file'])
        if filename.lower().endswith('.sdf'):
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
            
        elif (filename.lower().endswith('.smi') or 
                filename.lower().endswith('.ism') or 
                filename.lower().endswith('txt')):
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
        elif (filename.lower().endswith('.mol')):
            checks = request.POST.getlist('checks')            
            df = read_mol(request.FILES['file'], checks)            
            # Get table columns
            cols = list(df.columns)
            # Extract DataFrame rows
            rows = []
            for ix, row in df.iterrows():
                rows.append(zip(cols,list(row)))
            return render(request, 'table.html', {'PAGE_URL':PAGE_URL,
                                                  'rows':rows, 'cols':cols})
        elif (filename.lower().endswith('.mol2')):
            checks = request.POST.getlist('checks')            
            df = read_mol2(request.FILES['file'], checks)            
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
        
def api(request):
    """
    View that serves the API. SMILES have to be percent encoded
    """
    if request.method == 'GET':
        calcs = []        
        if 'calculations' in request.GET:
            if len(request.GET['calculations']) == 0:
                return JsonResponse({'calculations':['MW', 'logP', 'HBA', 'HBD', 
                                    'logS', 'PAINS']})
            elif 'smiles' in request.GET:
                calcs = request.GET['calculations'].split(',')
                if 'everything' in calcs:
                    calcs = ['MW', 'logP', 'HBA', 'HBD', 'logS', 'PAINS']
        if 'models' in request.GET and extra:
            if len(request.GET['models']) == 0:
                return JsonResponse(extra([], info=True))
            else:
                calcs.append('extra')
        
        if 'smiles' in request.GET and MolFromSmiles(request.GET['smiles']) is not None:
            df = read_smi_string(request.GET['smiles'], calcs)
            return JsonResponse(df.drop(['ROMol', '#'], axis=1).iloc[0].to_dict())     
        else:
            return HttpResponseBadRequest('Molecule could not be constructed from SMILES')
    else:
        return HttpResponseNotAllowed('Only GET')
