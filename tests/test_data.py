from yoke.annotations import RTKModules
from yoke.data.utils import BundledModules


def test_gbm_enum():
    """Test the GBM module enum is created correctly"""

    gbm = BundledModules.GBM
    assert gbm.name == "GBM"
    assert gbm.value.definition.exists()
    assert gbm.value.description.exists()
    assert gbm.value.hierarchy.exists()


def test_gbm_from_string():
    """Test that GBM modules can be created from string"""

    gbm = RTKModules.from_bundled("GBM")
    assert gbm.name == "GBM"
    assert gbm.definition.exists()
    assert gbm.description.exists()
    assert gbm.hierarchy.exists()
