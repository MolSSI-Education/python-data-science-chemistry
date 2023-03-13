from rdkit import Chem

class Dict2Mol():
    bond_types = {
        0: Chem.BondType.UNSPECIFIED,
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.QUADRUPLE,
        5: Chem.BondType.QUINTUPLE,
        6: Chem.BondType.HEXTUPLE,
        7: Chem.BondType.ONEANDAHALF,
        8: Chem.BondType.TWOANDAHALF,
        9: Chem.BondType.THREEANDAHALF,
        10: Chem.BondType.FOURANDAHALF,
        11: Chem.BondType.FIVEANDAHALF,
        12: Chem.BondType.AROMATIC,
        13: Chem.BondType.IONIC,
        14: Chem.BondType.HYDROGEN,
        15: Chem.BondType.THREECENTER,
        16: Chem.BondType.DATIVEONE,
        17: Chem.BondType.DATIVE,
        18: Chem.BondType.DATIVEL,
        19: Chem.BondType.DATIVER,
        20: Chem.BondType.OTHER,
        21: Chem.BondType.ZERO
    }
    mol_list = []
    names = []

    def __init__(self, **mol_set: dict):
        self.names = list(mol_set.keys())
        base_mol = Chem.Mol()
        for name in self.names:
            rw_mol = Chem.RWMol(base_mol)
            rw_mol.SetProp("_Name", name)
            mol = self.build_mol(rw_mol, **mol_set[name])
            self.mol_list += [mol]

    def build_mol(self, rw_mol, **mol_dict):
        for symbol in mol_dict["symbols"]:
            rw_mol.AddAtom(Chem.Atom(symbol))

        for atom1, atom2, bond in mol_dict["connectivity"]:
            bondtype = self.bond_types[bond]
            rw_mol.AddBond(atom1, atom2, bondtype)

        mol = rw_mol.GetMol()

        geom = Chem.Conformer(len(mol_dict["symbols"]))

        for index, coords in enumerate(mol_dict["geometry"]):
            geom.SetAtomPosition(
                index,
                coords
            )

        mol.AddConformer(geom)

        probs = None
        probs = Chem.DetectChemistryProblems(mol)
        if len(probs) == 0:
            pass
            s_flags = Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP|Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY|Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES|Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS|Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
            Chem.rdmolops.SanitizeMol(mol, sanitizeOps=s_flags)
        else:
            print("Chemistry problems detected in one or more molecules")
            print(probs)

        return mol
