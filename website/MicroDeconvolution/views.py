from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse

from .models import ImageUpload
from .forms import ImageForm, your_form_name, CustomRGB
from django.http import HttpResponse
from django.core.management import call_command
from . import RandomWalkScript

#def index(request):
#    return HttpResponse("Hello, world. You're at the MicroDeconvolution page.")

def index(request):
    return render(request, 'home.html')


def contact(request):
    return render(request, 'contact.html')

def tutorial(request):
    return render(request, 'tutorial.html')


def singleimage(request):
    # Handle file upload
    if request.method == 'POST':
        form = ImageForm(request.POST, request.FILES)
        if form.is_valid():
            newim = ImageUpload(imfile=request.FILES['imfile'])
            newim.save()
            #print(newim(type))

            file_path = 'images/uploaded/' + request.FILES['filename']
            #save_directory = '/home/griffin/Desktop/MicroDeconvolution/website/media/images/uploaded/'
            #RandomWalkScript.random_walk_segmentation(file_path, save_directory)


            # Redirect to the image list after POST
            return HttpResponseRedirect(reverse('singleimage'))
    else:
        form = ImageForm()  # A empty, unbound form

    # Load image for the list page
    images = ImageUpload.objects.all()

    # Render list page with the image and the form
    return render(
        request,
        'singleimage.html',
        {'images': images, 'form': form}
    )


def uploadpage(request):
    # Handle file upload
    if request.method == 'POST':
        form = ImageForm(request.POST, request.FILES)
        if form.is_valid():
            newim = ImageUpload(imfile=request.FILES['imfile'])
            newim.save()

            # Redirect to the image list after POST
            return HttpResponseRedirect(reverse('uploadpage'))
    else:
        form = ImageForm()  # A empty, unbound form

    # Load image for the list page
    images = ImageUpload.objects.all()

    # Render list page with the image and the form
    return render(
        request,
        'uploadpage.html',
        {'images': images, 'form': form}
    )


'''
import subprocess

def your_view_name(request):
    if request.method == 'GET':
        form = your_form_name()
        print('hello')
    else:
        if your_form_name.is_valid():
            print('hello')
            info = request.POST['info_name']
            output = script_function(info)
            #Here you are calling script_function,
            #passing the POST data for 'info' to it;
            return render(request, 'MicroDeconvolution/uploadpage.html', {
            'info': info,
            'output': output,
            })
        
    return render(request, 'MicroDeconvolution/uploadpage.html', {
        'form': form,
    })

def script_function( post_from_form ):
  print(post_from_form) #optional,check what the function received from the submit;
  return subprocess.check_call(['RandomWalkScript.py', post_from_form])
'''