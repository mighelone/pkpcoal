import pkp
import numpy as np

d = {'a': 23,
     'b': [0, 1, 2],
     'c': {
         'x': 1,
         'y': 1},
     'd': np.float(23),
     'e': np.arange(4)
     }


def test_clean_dict():
    clean_d = pkp.clean_dict(d)
    assert 'e' not in clean_d
    assert 'd' in clean_d


def test_clean_dict_list():
    clean_d = pkp.clean_dict(d, tolist=True)
    assert 'e' in clean_d
    assert isinstance(clean_d['d'], float)
    assert isinstance(clean_d['e'], list)


def test_clean_list():
    l = ['a', np.float(23.4)]
    clean_l = pkp.clean_list(l)
    assert isinstance(clean_l[1], float)

    clean_l = pkp.clean_list(tuple(l))
    assert isinstance(clean_l[1], float)
