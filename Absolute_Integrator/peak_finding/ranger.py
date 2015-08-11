import numpy as np
import scipy.signal
from scipy.ndimage.morphology import binary_erosion
from scipy.ndimage.morphology import white_tophat
from scipy.ndimage.filters import gaussian_filter

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

def normalise_dynamic_range(image):
    image -= image.min()
    image /= image.max()
    return image


    return image - gaussian_filter(image, filter_width)

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
    return 19

def get_end_search(image, end_search="auto"):
    im_dim = image.shape
    if end_search== "auto":
        return 2 * np.floor(( float(np.min(im_dim)) / 8) / 2) - 1
    else:
        return end_search

def fit_block(block, base_axis):
    A = np.vstack([base_axis**2 , base_axis , np.ones(base_axis.size)]).T
    h_profile = np.sum(block, axis=0)
    v_profile = np.sum(block, axis=1)
    solution_h = np.linalg.lstsq(A, np.log(h_profile))[0]
    solution_v = np.linalg.lstsq(A, np.log(v_profile))[0]

    y = -solution_v[1]/solution_v[0]/2.0
    x = -solution_h[1]/solution_h[0]/2.0
    height = ( h_profile.max() + v_profile.max() ) / 2.0
    spread = np.sqrt((np.abs(solution_h[0])+np.abs(solution_v[0])) / 4.0)

    return y, x, height, spread

# Feature identification section:
def filter_peaks(normalized_heights, spread, offset_radii, trial_size, sensitivity_threshold):

    # Normalise distances and heights:
    normalized_heights[normalized_heights < 0] = 0  # Forbid negative (concave) Gaussians.
    spread /= trial_size
    spread(spread > sqrt(2)) = sqrt(2) ;
    spread(spread == 0) = sqrt(2) ;
    offset_radii = offset_radii / trial_size
    offset_radii[offset_radii == 0] = 0.001  # Remove zeros values to prevent division error later.

    # Create search metric and screen impossible peaks:
    search_record = normalized_heights / offset_radii
    search_record /= 100.0
    search_record[search_record > 1] = 1
    search_record[spread < 0.5] = 0       # Invalidates negative Gaussian widths.
    search_record[spread > 1] = 0          # Invalidates Gaussian widths greater than a feature spacing.
    search_record[offset_radii > 1] = 0    # Invalidates Gaussian widths greater than a feature spacing.
    kernel = int(np.round(trial_size/3))
    if kernel % 2 == 0:
        kernel += 1
    search_record = scipy.signal.medfilt2d(search_record, kernel)  # Median filter to strip impossibly local false-positive features.
    search_record[search_record < sensitivity_threshold ] = 0   # Collapse improbable features to zero likelyhood.
    search_record[search_record >= sensitivity_threshold ] = 1  # Round likelyhood of genuine features to unity.
               
    # Erode regions of likely features down to points.
    search_record = binary_erosion(search_record, iterations=-1 )
    y, x = np.where(search_record==1)
    return np.vstack((y,x)).T  # Extract the locations of the identified features.



def peak_find(image,
              best_size="auto",
              refine_positions=False,
              sensitivity_threshold=33,
              start_search=3,
              end_search="auto",
              progress_object=None):
    """
    
    Parameters
    ----------
    refine_position : bool
        ddf
            
    """
    # TODO: best_size needs its auto-estimation routine
    trial_size = get_trial_size(best_size)

    # Removes slowly varying background from image to simplify Gaussian fitting.
    input_offset = white_tophat(image, 2*trial_size)

    # image dimension sizes, used for loop through image pixels
    m, n = get_data_shape(image)

    big = get_end_search(image, end_search)
            
    # Create blank arrays.
    heights        = np.empty(image.shape, dtype=np.float32)
    spreads         = np.empty(image.shape, dtype=np.float32)
    xs              = np.empty(image.shape, dtype=np.float32)
    ys              = np.empty(image.shape, dtype=np.float32)

        
    # Half of the trial size, equivalent to the border that will not be inspected.
    test_box_padding = int(( trial_size - 1 ) / 2.)

    # Coordinate set for X and Y fitting.  
    base_axis = np.arange(-test_box_padding, test_box_padding+1., dtype=np.float32)
    # Followed by the restoration progress bar:
    if progress_object is not None:
        progress_object.set_title("Identifying Image Peaks...")
        progress_object.set_position(0)
    for i in range(test_box_padding + 1 , m - ( test_box_padding + 1 )):
        currentStrip = input_offset[ i - test_box_padding : i + test_box_padding +1] 
        for j in range( test_box_padding + 1, n - ( test_box_padding + 1 )):
            I = currentStrip[:, j - test_box_padding : j + test_box_padding + 1]
            y, x, height, spread = fit_block(I, base_axis)
            ys[i, j] = y
            xs[i, j] = x
            heights[i, j] = height
            spreads[i, j] = spread
            
            if progress_object is not None:
                percentage_refined = (((trial_size-3.)/2.) / ((big-1.)/2.)) +  (((i-test_box_padding) / (m - 2*test_box_padding)) / (((big-1)/2)))  # Progress metric when using a looping peak-finding waitbar.
                progress_object.set_position(percentage_refined)
    # normalize peak heights
    heights = heights / ( np.max(input_offset) - np.min(input_offset) ) 
    # normalize fitted Gaussian widths
    spreads = spreads / trial_size
    offset_radii = np.sqrt(ys**2 + xs**2)  # Calculate offset radii.
    return filter_peaks(heights, spreads, offset_radii, trial_size, sensitivity_threshold)
