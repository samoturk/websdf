from django.shortcuts import render

def home(request):
    '''
    Define home page view. It is loaded by urls.py.
    '''
    
    return render(request, 'index.html')