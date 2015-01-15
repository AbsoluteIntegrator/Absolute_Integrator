import numpy as np
import skimage.feature
import skimage.filter

from Absolute_Integrator.pre_processing import background_subtraction
from Absolute_Integrator.pre_processing import peak_width_estimation


# dictionary describing options available to tune this algorithm
options = {
    "peak_width": {"purpose": "The estimate of the peak size, in pixels.  If 'auto', attempts to determine automatically.  Otherwise, this should be an integer.",
                   "default": "auto"},
    "blur": {"purpose": "optionally blur image before attempting to find peaks",
             "default": True},
    "subtract_background": {"purpose": "optionally subtract background around peaks to improve location",
                            "default": True},
}

def peak_find(image, peak_width="auto", blur=True, subtract_background=True):
    if blur:
        # first blur the image so we don't find spurious peaks
        image = skimage.filter.gaussian_filter(image, 3)

    if subtract_background:
        # subtracting the background makes the peaks much easier to locate
        image = background_subtraction.subtract_background(image)

    if peak_width == "auto":
        # estimate the peak width to know what size peaks to look for
        peak_width = peak_width_estimation.estimate_peak_width(image)

    # people think of peaks in terms of diameter, but sigma is expressed in radius
    peak_width /= (2 * np.sqrt(2.0))

    # search in a neighborhood around the estimated peak size
    min_peak_size = peak_width - 0.25 * peak_width
    max_peak_size = peak_width + 0.25 * peak_width
    peaks = skimage.feature.blob_log(image, min_sigma=min_peak_size, max_sigma=max_peak_size, num_sigma=4)

    # make third column of data be approximate peak diameter
    peaks[:, 2] = np.array(peaks[:, 2] * 2 * np.sqrt(2.0), dtype=np.int64)
    return peaks