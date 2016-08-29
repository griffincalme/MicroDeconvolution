from django import forms


class ImageForm(forms.Form):
    imfile = forms.ImageField(
        label='Select a file',

    )

class your_form_name(forms.Form):
    error_css_class = 'error' #custom css for form errors - ".error";
    required_css_class = 'required' #custom css for required fields - ".required";
    info_text = forms.CharField()

#class StainRGB(forms.Form):


class CustomRGB(forms.Form):
    red_val = forms.ChoiceField(choices=[(x, x) for x in range(0,255)])
    green_val = forms.ChoiceField(choices=[(x, x) for x in range(0,255)])
    blue_val = forms.ChoiceField(choices=[(x, x) for x in range(0,255)])
    label = forms.CharField(widget=forms.TextInput(attrs={'placeholder': 'Required', 'size': 20}), error_messages={'required': 'Required'})