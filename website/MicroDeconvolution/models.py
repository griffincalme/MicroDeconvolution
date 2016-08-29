from django.db import models


class ImageUpload(models.Model):
    imfile = models.ImageField(upload_to='images/uploaded/', default='NULL')


'''
class ProcessParameters(models.Model):
    STAINS = (
        ('DAB', ' a'),
        ('Hematoxylin', 'b '),
    )
    stain = models.char
'''