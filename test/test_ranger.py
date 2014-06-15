__author__ = 'msarahan'

import pytest
import numpy as np

import pattern
import Absolute_Integrator.peak_finding.ranger as ranger


def test_auto_end_search():
    image = np.zeros((256, 256))
    end_search = ranger.get_end_search(image)
    assert end_search == 2 * np.floor((np.min(image.shape) / 8.0) / 2) - 1


def test_auto_trial_size():
    offset = 32
    image = pattern.get_test_pattern((256, 256), offset=offset)
    trial_size = ranger.get_trial_size(image)
    # assert that we're within 10% of the known spacing between peaks
    assert 0.9*offset < trial_size < 1.1 * offset


def test_image_size():
    image = np.zeros((512, 256))
    shape = ranger.get_data_shape(image)
    # our data shape is the flip of numpy order.  We should straighten this out eventually, but test for it here.
    assert shape == image.shape[::-1]

def test_estimate_gaussian_parameters():
    data_size = 15
    height = 232
    # slightly off-center
    center = 1.0
    sigma = 4.2
    # generate a 1D gaussian of known parameters:
    gaussian = lambda x: height*np.exp(-(((center-x)/sigma)**2))
    base_axis = np.arange(-data_size/2, data_size/2)
    gaussian = gaussian(base_axis)
    # assert that we've measured it properly
    measured_center, measured_sigma, measured_height = ranger.estimate_1D_Gaussian_parameters(gaussian, base_axis)
    # center within half a pixel
    assert np.isclose(measured_center, center, 0.5)
    # sigma within half a pixel
    assert np.isclose(measured_sigma, sigma, 0.5)
    # height within 10%
    assert np.isclose(measured_height, height, .1*height)


