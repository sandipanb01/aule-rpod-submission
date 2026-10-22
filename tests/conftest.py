import pytest

from aule_rpod.config import load_config
from aule_rpod.scenario import Scenario


@pytest.fixture(scope="session")
def cfg():
    return load_config()


@pytest.fixture(scope="session")
def sc(cfg):
    return Scenario(cfg)
