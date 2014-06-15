import numpy as np
import numpy.linalg
import scipy.signal
from scipy.ndimage.morphology import binary_erosion

# dictionary describing options available to tune this algorithm
options = {
    "best_size": {"purpose": "Estimate of the distance between peaks, in pixels.  If 'auto', attempts to determine automatically.  Otherwise, this should be an odd integer.",
                  "default": "auto",
                  "type": "int",
                  "has_auto": True},
    "refine_positions": {"purpose": "TODO",
                         "default": False,
                         "type": "bool"},
    "sensitivity_threshold": {"purpose": "TODO",
                              "default": 0.34,
                              "type": "float"},
    "start_search": {"purpose": "TODO",
                     "default": 3,
                     "type": "int"},
    "end_search": {"purpose": "TODO",
                   "default": "auto",
                   "type": "int",
                   "has_auto": True},
    "progress_object": {"purpose": "Object used to present a progress bar to the user.  For definition, see UI_interface folder.",
                        "default": None},
}

def estimate_1D_Gaussian_parameters(data, axis):
    center = np.sum(axis * data) / np.sum(data)
    sigma = np.sqrt(np.abs(np.sum((axis - center) ** 2 * data) / np.sum(data)))
    height = data.max()
    return center, sigma, height

def remove_dc_offset(image):
    """ Removes 'DC offset' from image to simplify Gaussian fitting. """
    return (image - image.min()).astype("float32")

def get_data_shape(image):
    """ Returns data shape as (columns, rows).  
    Note that this is opposite of standard Numpy/Hyperspy notation.  Presumably,
    this is to help Lewys keep X and Y in order because Matlab uses column-major indexing.
    """
    n, m = image.shape
    return m, n

def get_trial_size(image, best_size="auto"):
    """ TODO: automatically estimate best box size """
    return 10

def get_end_search(image, end_search="auto"):
    im_dim = image.shape
    # Search should terminate around 1/8th of the smallest image dimension. Any larger has no physcial meaning
    # (fewer than 4 periodic features at the Nyquist limit).  Must also be rounded to the neareast ODD number.    
    if end_search == "auto":
        return 2 * np.floor((float(np.min(im_dim)) / 8) / 2) - 1
    else:
        return end_search

def fit_block(block, base_axis):
    x, sx, hx = estimate_1D_Gaussian_parameters(np.sum(block, axis=0), base_axis) # The horizontal offset refinement.
    y, sy, hy = estimate_1D_Gaussian_parameters(np.sum(block, axis=1), base_axis) # The vertical offset refinement.
                
    # use base_axis length as way of not passing trial_size as parameter
    height = (hx+hy) /(2*len(base_axis))  # Calculates the height of the fitted Gaussian.
    spread = 2.3548 * np.sqrt(sx ** 2 + sy ** 2)  # 2D FWHM
    return y, x, height, spread

# Feature identification section:
def filter_peaks(normalized_heights, spread, offset_radius, trial_size, sensitivity_threshold):
    normalized_heights[normalized_heights < 0] = 0  # Forbid negative (concave) Gaussians.
    offset_radius[offset_radius == 0] = np.nan  # Remove zeros values to prevent division error later.
    # Create search metric and screen impossible peaks:
    search_record = normalized_heights / offset_radius
    search_record[search_record > 1] = 1 
    search_record[search_record < 0] = 0 
    search_record[spread < 0.05] = 0       # Invalidates negative Gaussian widths.
    search_record[spread > 1] = 0          # Invalidates Gaussian widths greater than a feature spacing.
    search_record[offset_radius > 1] = 0    # Invalidates Gaussian widths greater than a feature spacing.
    kernel = int(np.round(trial_size/3))
    if kernel % 2 == 0:
        kernel += 1
    search_record = scipy.signal.medfilt2d(search_record, kernel)  # Median filter to strip impossibly local false-positive features.
    search_record[search_record < sensitivity_threshold] = 0   # Collapse improbable features to zero likelyhood.
    search_record[search_record >= sensitivity_threshold] = 1  # Round likelyhood of genuine features to unity.
               
    # Erode regions of likely features down to points.
    search_record = binary_erosion(search_record, iterations=-1)
    y, x = np.where(search_record == 1)
    return np.vstack((y, x)).T  # Extract the locations of the identified features.



def run(image,
        best_size="auto",
        refine_positions=False,
        sensitivity_threshold=0.34,
        start_search=3,
        end_search="auto",
        progress_object=None):
    """
    
    Parameters
    ----------
    refine_position : bool
        ddf
            
    """
    # Removes 'DC offset' from image to simplify Gaussian fitting.
    input_offset = remove_dc_offset(image)

    # image dimension sizes, used for loop through image pixels
    m, n = get_data_shape(image)

    big = get_end_search(image, end_search)
            
    # TODO: best_size needs its auto-estimation routine
    trial_size = get_trial_size(image, best_size)

    # Create blank arrays.
    heights        = np.empty(image.shape) 
    spread         = np.empty(image.shape)
    x              = np.empty(image.shape)
    y              = np.empty(image.shape)

    # Half of the trial size, equivalent to the border that will not be inspected.
    test_box_padding = int((trial_size - 1) / 2.)

    # Coordinate set for X and Y fitting.  
    base_axis = np.arange(-test_box_padding, test_box_padding+1)
    # Followed by the restoration progress bar:
    if progress_object is not None:
        progress_object.set_title("Identifying Image Peaks...")
        progress_object.set_position(0)
    for i in xrange(test_box_padding + 1, m - (test_box_padding + 1)):
        currentStrip = input_offset[i - test_box_padding: i + test_box_padding + 1]
        for j in xrange(test_box_padding + 1, n - (test_box_padding + 1 )):
            I = currentStrip[:, j - test_box_padding: j + test_box_padding + 1]
            y[i, j], x[i, j], heights[i, j], spread[i, j] = fit_block(I, base_axis)
            
            if progress_object is not None:
                percentage_refined = (((trial_size-3.)/2.) / ((big-1.)/2.)) + (((i-test_box_padding) / (m - 2 * test_box_padding)) / (((big - 1) / 2)))  # Progress metric when using a looping peak-finding waitbar.
                progress_object.set_position(percentage_refined)
    # normalize peak heights
    heights = heights / (np.max(input_offset) - np.min(input_offset))
    # normalize fitted Gaussian widths
    spread /= trial_size
    offset_radius = np.sqrt(y**2 + x**2)  # Calculate offset radii.
    offset_radius /= trial_size
    return filter_peaks(heights, spread, offset_radius, trial_size, sensitivity_threshold) 
