def split_pdb_by_chain(pdb_str):
    """Split a PDB string into a list of PDB strings, one for each chain."""
    pdb_list = []
    chain = None
    chain_str = ""
    for line in pdb_str.splitlines():
        if line.startswith("ATOM"):
            if chain is None:
                chain = line[21]
            elif line[21] != chain:
                pdb_list.append(chain_str)
                chain_str = ""
                chain = line[21]
            chain_str += line + "\n"
        else:
            continue
    pdb_list.append(chain_str)
    return pdb_list
