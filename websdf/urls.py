from django.conf.urls import url
#from django.contrib import admin
from websdf.settings import APP_NAME
from websdf import views

urlpatterns = [
    url('^' + APP_NAME + r'$', views.home, name='home'),
    url('^' + APP_NAME + r'sdf', views.upload_file, name='upload_file'),
    url('^' + APP_NAME + r'api', views.api, name='api'),
]
