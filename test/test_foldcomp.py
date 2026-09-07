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


def test_decompress_batch_cpu_handles_raw_container_fragments(pytestconfig):
    pdb = pytestconfig.rootpath.joinpath("test/test.pdb").read_bytes()
    raw_container = foldcomp.compress("raw", pdb, max_backbone_rmsd=0.0)
    assert raw_container.startswith(b"FCZC")
    [(name, text)] = foldcomp.decompress_batch([raw_container], use_gpu=False)
    assert name == "raw"
    assert any(line.startswith("ATOM") for line in text.splitlines())
