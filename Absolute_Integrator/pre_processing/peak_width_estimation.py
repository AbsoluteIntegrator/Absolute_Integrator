import numpy as np
from copy import deepcopy


def find_bounds(one_dim_array):
    """Attempt to find peak bounds given a 1D strip containing a peak at its center"""
    zero_crossings = np.where(np.diff(np.sign(np.diff(one_dim_array))))[0]
    # reuse the cx variable to represent the middle of our sub-window.
    # they're 0-based indexes, so adjust for that...
    cx = one_dim_array.shape[0] / 2 - 1
    if len(zero_crossings) < 2 or min(zero_crossings) > cx or max(zero_crossings) < cx:
            return 0
    # find the zero crossing closest to the left of the peak.
    left = zero_crossings[zero_crossings < cx-1][-1]
    # First, chuck any values to the left of our peak.
    peak_right_crossings = zero_crossings[zero_crossings > cx+1]
    # The leftmost value is the closest one.
    if peak_right_crossings.size > 0:
        right = peak_right_crossings[0]
    else:
        right = one_dim_array.size
    return left, right


def estimate_peak_width(image, window_size=64):
    """Given input image, locate highest peak and estimate its width.

    Width is estimated by examining derivatives of strips vertically and horizontally from the peak.

    window_size is the window around the highest peak used to estimate the width.
    """
    # make a copy of the image so that we don't alter people's data
    tmp_image = deepcopy(image)
    # apply a quick mask so that we don't end up on the edges.
    tmp_image[0:window_size / 2] = 0
    tmp_image[-window_size/2:] = 0
    tmp_image[:, 0:window_size/2] = 0
    tmp_image[:, -window_size/2:] = 0
    # find the highest point in the image.
    k = np.argmax(tmp_image)
    cx, cy = np.unravel_index(k, image.shape)
    # Pick out the middle row through that highest point. Make sure datatype
    # is signed so we can have negative numbers
    window_offset = window_size / 2
    tmp_row = np.array(tmp_image[cy, cx-window_offset:cx+window_offset], dtype=np.float64)
    tmp_col = np.array(tmp_image[cy-window_offset:cy+window_offset, cx], dtype=np.float64)
    # Detect where our derivative switches sign
    left, right = find_bounds(tmp_row)
    top, bottom = find_bounds(tmp_col)
    # width is average of vertical and horizontal widths
    width = np.mean([right-left, bottom-top])

    return width