from pyket.parse_input import read_photon_bath

def test_parse_input():
    with open("fake_input.txt") as f:
        results = read_photon_bath(f,True)
        assert len(results) == 2
        assert results['5'] == '6'
        assert results['1.2'] == '2.5'
