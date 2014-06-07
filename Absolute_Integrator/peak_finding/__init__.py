import ranger
import template_match

methods={"ranger": ranger,
         "template_match": template_match,
         }

# default peak finding method is ranger.  How should people be able to specify the method?
default = "Ranger"

def peak_find(image, method=default, **options):
    # look up which method to use from the dict of methods
    if method.lower() in methods:
        method = methods[method.lower()]
    else:
        raise ValueError("Peak finding method {:s} not recognized.  Available methods: {:s}".format(
            method, str(methods.keys)))
    return method.peak_find(image, **options)
