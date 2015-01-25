from django.conf.urls import patterns, include, url
#from django.contrib import admin

urlpatterns = patterns('',
     url(r'^$', 'websdf.views.home', name='home'),
     url(r'^sdf', 'websdf.views.upload_file', name='upload_file'),
)
