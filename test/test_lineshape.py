import pytest
from pyket.lineshapes import lorentzian

def test_lorentzian_max():
    p0 = 2.3
    width = 3.0
    actual_max = lorentzian(p0,p0,width)
    expected_max = pytest.approx(1.0,abs=0.01)
    assert expected_max == actual_max

def test_lorentzian_symmetric():
    p0 = -5.6
    width = .66
    delta_x = 20
    left_val = lorentzian(-delta_x+p0, p0, width)
    expected_right_val = pytest.approx(left_val,abs=0.0001)
    actual_right_val = lorentzian(delta_x+p0,p0,width)
    assert expected_right_val == actual_right_val
