import numpy as np
import numpy.linalg
import scipy.signal
from scipy.ndimage.morphology import binary_erosion

# dictionary describing options available to tune this algorithm
options = {
    "best_size":{"purpose":"The estimate of the peak size, in pixels.  If 'auto', attempts to determine automatically.  Otherwise, this should be an integer.",
                 "default":"auto"},
    "refine_positions":{"purpose":"TODO",
                        "default":False},
    "sensitivity_threshold":{"purpose":"TODO",
                             "default":0.34},
    "start_search":{"purpose":"TODO",
                    "default":3},
    "end_search":{"purpose":"TODO",
                  "default":"auto"},
    "progress_object":{"purpose":"Object used to present a progress bar to the user.  For definition, see UI_interface folder.",
                       "default":None},
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
    im_dim = image.shape[::-1]
    m, n = im_dim
    return m, n

def get_trial_size(image, best_size="auto"):
    """ TODO: automatically estimate best box size """
    return 10

def get_end_search(image, end_search="auto"):
    im_dim = image.shape
    if end_search== "auto":
        return 2 * np.floor(( float(np.min(im_dim)) / 8) / 2) - 1
    else:
        return end_search

def fit_block(block, base_axis):
    x, sx, hx = estimate_1D_Gaussian_parameters(np.sum(block, axis=0), base_axis) # The horizontal offset refinement.
    y, sy, hy = estimate_1D_Gaussian_parameters(np.sum(block, axis=1), base_axis) # The vertical offset refinement.
                
    # use base_axis length as way of not passing trial_size as parameter
    height = (hx+hy) /(2*len(base_axis))  # Calculates the height of the fitted Gaussian.
    spread =   2.3548 * np.sqrt(sx ** 2 + sy ** 2)  # 2D FWHM
    return height, spread

# Feature identification section:
def filter_peaks(peaks, spread, sensitivity_threshold):
    normalisedPeak = peaks / ( np.max(inputOffset) - np.min(inputOffset) )  # Make use of peak knowledge:
    normalisedPeak[normalisedPeak < 0] = 0  # Forbid negative (concave) Gaussians.
    spread = spread / trialSize          # Normalise fitted Gaussian widths.
    offsetRadius = np.sqrt( (horz_offset)**2 + (vert_offset)**2 )  # Calculate offset radii.
    offsetRadius = offsetRadius / trialSize 
    offsetRadius[offsetRadius == 0] = 0.001  # Remove zeros values to prevent division error later.
    # Create search metric and screen impossible peaks:
    search_record = normalisedPeak / offsetRadius
    search_record[search_record > 1] = 1 
    search_record[search_record < 0] = 0 
    search_record[spread < 0.05] = 0       # Invalidates negative Gaussian widths.
    search_record[spread > 1] = 0          # Invalidates Gaussian widths greater than a feature spacing.
    search_record[offsetRadius > 1] = 0    # Invalidates Gaussian widths greater than a feature spacing.
    kernel = int(np.round(trialSize/3))
    if kernel % 2 == 0:
        kernel += 1
    search_record = scipy.signal.medfilt2d(search_record, kernel)  # Median filter to strip impossibly local false-positive features.
    sensitivityThreshold = 0.34            # This is an Admin tunable parameter that is defined here within the core file.
    search_record[search_record < sensitivity_threshold ] = 0   # Collapse improbable features to zero likelyhood.
    search_record[search_record >= sensitivity_threshold ] = 1  # Round likelyhood of genuine features to unity.
               
    # Erode regions of likely features down to points.
    search_record = binary_erosion(search_record, iterations=-1 )     
    # [point_coordinates(:,2),point_coordinates(:,1)] = np.where(search_record == 1)  # Extract the locations of the identified features.
    #    return search_record


def peak_find(image,
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
    inputOffset = remove_dc_offset(image)

    # image dimension sizes, used for loop through image pixels
    m, n = get_data_shape(image)

    big = get_end_search(end_search)
            
    # TODO: best_size needs its auto-estimation routine
    trialSize = get_trial_size(best_size)

    # Create blank arrays.
    peak        = np.empty(image.shape) 
    spread      = np.empty(image.shape)
        
    # Half of the trial size, equivalent to the border that will not be inspected.
    test_box_padding = int(( trialSize - 1 ) / 2.)

    # Coordinate set for X and Y fitting.  
    base_axis = np.arange(-test_box_padding, test_box_padding)
    # Followed by the restoration progress bar:
    if progress_object is not None:
        progress_object.set_title("Identifying Image Peaks...")
        progress_object.set_position(0)
    for i in xrange(test_box_padding + 1 , m - ( test_box_padding + 1 )):
        currentStrip = inputOffset[ i - test_box_padding : i + test_box_padding] 
        for j in xrange( test_box_padding + 1, n - ( test_box_padding + 1 )):
            I = currentStrip[:, j - test_box_padding : j + test_box_padding]
            peak[i,j], spread[i,j] = fit_block(I, base_axis)
            
            if progress_object is not None:
                percentage_refined = ( ((trialSize-3.)/2.) / ((big-1.)/2.) ) +   ( ( (i-test_box_padding) / (m - 2*test_box_padding) ) / (((big-1)/2)))  # Progress metric when using a looping peak-finding waitbar.
                progress_object.set_position(percentage_refined)
    return filter_peaks(peak, spread, sensitivityThreshold) 
