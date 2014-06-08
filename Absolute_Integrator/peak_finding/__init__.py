import ranger
import template_match

methods={"ranger": ranger,
         "template_match": template_match,
         }

# default peak finding method is ranger.
default_method = "Ranger"

def list_methods():
    print methods.keys

def list_options(method):
    """
    List the available options for a given peak finding routine.

    Parameters:
    -----------
    method : string
         The string identifier of the method to show options for.  To show
         all available methods, call the list_methods function.

    Returns:
    dict : keys are options, values are descriptions of those options.
    """
    if method.lower() in methods:
        return methods[method.lower()].options

def peak_find(image, method=default_method, **options):
    """
    Executes a peak finding method on given data.

    Parameters:
    -----------
    image : ndarray
        The input image to find peaks on
    method : string
        The string identifier of the method to use.  To show what
        methods are available, call the list_methods function.
    **options : arbitrary key-specified options to be passed to peak finder.
        To find what options are available, call the list_options function,
        with the string identifier of the method you'd like to use.

    Returns:
    --------
    ndarray : n x 2, with n being the number of peaks found.  The positions
        are in numpy order (Y, X), so that they can be fed directly to 
        coordinate extraction routines.
    """
    # look up which method to use from the dict of methods
    if method.lower() in methods:
        method = methods[method.lower()]
    else:
        raise ValueError("Peak finding method {:s} not recognized.  Available methods: {:s}".format(
            method, str(methods.keys)))
    return method.peak_find(image, **options)