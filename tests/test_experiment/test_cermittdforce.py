from mrfmsim.experiment import CermitTDGroup_force
import pytest


"""Test CermitTDGroup_force module in mrfmsim.experiment

There is no need to test the structure of the simulator, so we just want to import it to make sure it is working properly. The tests for the simulator are in the test_experiment/test_cermittd.py file.

"""

CermitTDForce_mul = CermitTDGroup_force.experiments["CermitTD_force_multipulse"]
CermitTDForce_one = CermitTDGroup_force.experiments["CermitTD_force_onepulse"]


def test_cermittdforce_multipulse_import():
    """Test that the multipulse experiment can be successfully referenced."""

    assert CermitTDForce_mul is not None


def test_cermittdforce_onepulse_import():
    """Test that the onepulse experiment can be successfully referenced."""

    assert CermitTDForce_one is not None
