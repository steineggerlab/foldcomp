import foldcomp
import pytest
import shutil
from pathlib import Path


def test_decompress(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/test.pdb"), "rb") as f:
        data = f.read().decode("utf-8")
        print(foldcomp.decompress(foldcomp.compress("test", data)))


def test_decompress_mmcif(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/test.pdb"), "rb") as f:
        data = f.read().decode("utf-8")
    fcz = foldcomp.compress("test", data)
    name, mmcif = foldcomp.decompress(fcz, format="mmcif")
    assert isinstance(name, str)
    assert isinstance(mmcif, str)
    assert mmcif.startswith("data_")
    assert "_atom_site.Cartn_x" in mmcif


def test_decompress_invalid_format(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/test_af.fcz"), "rb") as f:
        fcz = f.read()
    with pytest.raises(ValueError, match="format must be one of"):
        foldcomp.decompress(fcz, format="xyz")


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


def test_compress_multichain_requires_split(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/multichain.pdb"), "rb") as f:
        data = f.read().decode("utf-8")
    with pytest.raises(foldcomp.error, match="Multiple chains found"):
        foldcomp.compress("multichain.pdb", data)


def test_compress_split_multichain(pytestconfig):
    with open(pytestconfig.rootpath.joinpath("test/multichain.pdb"), "rb") as f:
        data = f.read().decode("utf-8")
    chunks = foldcomp.compress("multichain.pdb", data, split=True)
    assert isinstance(chunks, list)
    assert len(chunks) > 1
    first_name, first_fcz = chunks[0]
    assert first_name.endswith(".fcz")
    assert isinstance(first_fcz, bytes)
    name, pdb = foldcomp.decompress(first_fcz)
    assert isinstance(name, str)
    assert isinstance(pdb, str)
    assert "ATOM" in pdb


def _copy_example_db(rootpath, tmp_path):
    source = rootpath.joinpath("test/example_db")
    target = tmp_path.joinpath("example_db_copy")
    for suffix in ["", ".dbtype", ".index", ".lookup", ".source"]:
        src = Path(str(source) + suffix)
        if src.exists():
            shutil.copy(src, Path(str(target) + suffix))
    return target


def test_open_merge_fragments_lookup_grouping(pytestconfig, tmp_path):
    dbpath = _copy_example_db(pytestconfig.rootpath, tmp_path)
    lookup_path = Path(str(dbpath) + ".lookup")
    lines = lookup_path.read_text().splitlines()
    first_cols = lines[0].split("\t")
    second_cols = lines[1].split("\t")
    first_cols[1] = "merged_entry"
    second_cols[1] = "merged_entry"
    lines[0] = "\t".join(first_cols)
    lines[1] = "\t".join(second_cols)
    lookup_path.write_text("\n".join(lines) + "\n")

    with foldcomp.open(str(dbpath)) as raw_db:
        raw_len = len(raw_db)

    with foldcomp.open(str(dbpath), merge_fragments=True) as merged_db:
        assert len(merged_db) == raw_len - 1
        found = False
        for i in range(len(merged_db)):
            name, pdb = merged_db[i]
            if name == "merged_entry":
                assert isinstance(pdb, str)
                assert "ATOM" in pdb
                source = merged_db.source_indices(i)
                assert len(source) == 2
                found = True
                break
        assert found

    with foldcomp.open(
        str(dbpath), ids=["merged_entry"], merge_fragments=True
    ) as merged_filtered:
        assert len(merged_filtered) == 1
        name, pdb = merged_filtered[0]
        assert name == "merged_entry"
        assert "ATOM" in pdb
        assert len(merged_filtered.source_indices(0)) == 2

    with foldcomp.open(
        str(dbpath), ids=["merged_entry"], merge_fragments=True, format="mmcif"
    ) as merged_mmcif:
        assert len(merged_mmcif) == 1
        name, mmcif = merged_mmcif[0]
        assert name == "merged_entry"
        assert isinstance(mmcif, str)
        assert mmcif.startswith("data_")
        assert mmcif.count("data_") == 1
        assert "_atom_site.Cartn_x" in mmcif


def test_open_merge_fragments_requires_decompress(pytestconfig, tmp_path):
    dbpath = _copy_example_db(pytestconfig.rootpath, tmp_path)
    with pytest.raises(TypeError, match="merge_fragments requires decompress=True"):
        foldcomp.open(str(dbpath), merge_fragments=True, decompress=False)


def test_source_indices_requires_merge(pytestconfig):
    with foldcomp.open(str(pytestconfig.rootpath.joinpath("test/example_db"))) as db:
        with pytest.raises(TypeError, match="source_indices is available only"):
            db.source_indices(0)


def test_open_invalid_format(pytestconfig):
    with pytest.raises(ValueError, match="format must be one of"):
        foldcomp.open(
            str(pytestconfig.rootpath.joinpath("test/example_db")), format="xyz"
        )
