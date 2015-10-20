from django.conf.urls import patterns, include, url
#from django.contrib import admin
from websdf.settings import APP_NAME

urlpatterns = patterns('',
    url('^' + APP_NAME + r'$', 'websdf.views.home', name='home'),
    url('^' + APP_NAME + r'sdf', 'websdf.views.upload_file', name='upload_file'),
)
