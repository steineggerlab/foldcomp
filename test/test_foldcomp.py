import foldcomp
from pathlib import Path


def test_decompress(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/test.pdb"), "rb") as f:
        data = f.read().decode("utf-8")
        print(foldcomp.decompress(foldcomp.compress("test", data)))


def test_open_db_all(pytestconfig):
    path = Path(pytestconfig.rootpath.joinpath("test/example_db"))
    with foldcomp.open(path) as db:
        for i in db:
            print(i)


def test_open_db_ids(pytestconfig):
    path = Path(pytestconfig.rootpath.joinpath("test/example_db"))
    with foldcomp.open(path, ids=["d1asha_", "d1it2a_"]) as db:
        for i in db:
            print(i)


def test_open_db_str(pytestconfig):
    with foldcomp.open(str(pytestconfig.rootpath.joinpath("test/example_db"))) as db:
        pass
