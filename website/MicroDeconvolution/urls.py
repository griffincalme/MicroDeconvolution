from django.conf.urls import url
from . import views

urlpatterns = [
    url(r'^$', views.index, name='index'),
    url(r'^contact/', views.contact, name='contact'),
    url(r'^tutorial/', views.tutorial, name='tutorial'),
    url(r'^singleimage/$', views.singleimage, name='singleimage'),
    url(r'^uploadpage/$', views.uploadpage, name='uploadpage')
]