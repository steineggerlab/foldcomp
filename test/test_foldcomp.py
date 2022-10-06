import foldcomp
from pathlib import Path


def test_decompress(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/test.pdb"), "rb") as f:
        data = f.read().decode("utf-8")
        print(foldcomp.decompress(foldcomp.compress("test", data)))


def test_open_db(pytestconfig):
    path = Path(pytestconfig.rootpath.joinpath("test/example_db"))
    ids = ["d1asha_", "d1it2a_"]
    with foldcomp.open(path) as db:
        for i in db:
            print(i)


def test_open_db_str(pytestconfig):
    ids = ["d1asha_", "d1it2a_"]
    with foldcomp.open(str(pytestconfig.rootpath.joinpath("test/example_db"))) as db:
        pass
