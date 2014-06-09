"""
Tests for `peak_finding` module.
"""
import pytest
import numpy as np

import Absolute_Integrator.peak_finding as peak_finding
import pattern

def test_list_methods():
    assert len(peak_finding.list_methods())>0
        
def test_list_options():
    methods = peak_finding.list_methods()
    # if any methods are missing the options object, this will raise an exception.
    for method in methods:
        peak_finding.list_options(method)

def test_methods():
    data, npeaks = pattern.get_test_pattern((256, 256))
    methods = peak_finding.list_methods()
    for method in methods:
        peaks = peak_finding.peak_find(data, method=method)
        assert peaks.shape==(npeaks,2)

def test_invalid_method():
    with pytest.raises(ValueError):
        peak_finding.peak_find(np.zeros((10,10)), method="Invalid")
