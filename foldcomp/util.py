def split_pdb_by_chain(pdb_str):
    """Split a PDB string into a list of PDB strings, one for each chain."""
    if isinstance(pdb_str, (bytes, bytearray)):
        pdb_str = pdb_str.decode("utf-8")
    pdb_list = []
    chain = None
    chain_str = ""
    for line in pdb_str.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            line_chain = line[21] if len(line) > 21 else ""
            if chain is None:
                chain = line_chain
            elif line_chain != chain:
                pdb_list.append(chain_str)
                chain_str = ""
                chain = line_chain
            chain_str += line + "\n"
        else:
            continue
    if chain is not None:
        pdb_list.append(chain_str)
    return pdb_list
