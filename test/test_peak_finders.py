"""
Tests for `Absolute_Integrator` module.
"""
import pytest
from Absolute_Integrator import peak_finding
import pattern

class TestPeakFinding(object):
    data=None

    @classmethod
    def setup_class(cls):
        self.data = pattern.get_test_pattern((256, 256))

    def test_ranger(self):
        peaks = peak_finding.peak_find(self.data, method="Ranger")
        assert peaks.shape==(20,2)

    @classmethod
    def teardown_class(cls):
        pass
