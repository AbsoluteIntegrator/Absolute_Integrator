# next few lines automatically collect modules for whatever methods you've implemented.
__all__ = ['ranger', 'template_match']
# Don't modify the line above, or this line!
import automodinit
automodinit.automodinit(__name__, __file__, globals())
del automodinit

# OK to modify from here down.
# create a dict of our modules to show what options are available to people
__methods = dict()
for module in __all__:
    __methods[module] = eval(module)
del module

# default peak finding method is ranger.
default_method = "Ranger"

def list_methods():
    return __methods.keys()

def list_options(method):
    """
    List the available options for a given peak finding routine.

    Parameters:
    -----------
    method : string
         The string identifier of the method to show options for.  Case insensitive. To show
         all available methods, call the list_methods function.

    Returns:
    dict : keys are options, values are descriptions of those options.
    """
    method = method.lower()
    if method in __methods:
        return __methods[method].options
    else:
        raise ValueError("Peak finding method {:s} not recognized.  Available methods: {:s}".format(
            method, str(list_methods())))

def peak_find(image, method=default_method, **options):
    """
    Executes a peak finding method on given data.

    Parameters:
    -----------
    image : ndarray
        The input image to find peaks on
    method : string
        The string identifier of the method to use.  Case insensitive.  To show what
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
    method = method.lower()
    if method in __methods:
        method = __methods[method]
    else:
        raise ValueError("Peak finding method {:s} not recognized.  Available methods: {:s}".format(
            method, str(list_methods())))
    return method.peak_find(image, **options)
