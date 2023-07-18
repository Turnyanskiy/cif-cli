import os
import pytest
from cif_cli.transformations import transform
import Bio


@pytest.mark.parametrize('filepath, expected', [
    ('./test_transformation/7g93.cif', True),
    ('./test_transformation/7g93.pdb', True),
    ('./test_transformation/7g93', False),
])
def test_get_model(filepath, expected):
    assert isinstance(transform.get_model(filepath), Bio.PDB.Model.Model) == expected


def test_save_model():
    filepath = './test_transformation/7g93.cif'
    model = transform.get_model(filepath)
    transform.save_model(model, filepath)
    os.system('rm -r 7g93_transformed.cif')
