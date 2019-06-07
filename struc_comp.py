# /usr/bin/env/ python3
## Find Closest residues to phosphate atoms
import Bio.PDB
import Bio
from Bio.PDB import PDBIO
import math

upath = "/home/jonah/Dropbox (ASU)/"
# macpath = "/Users/Amber/Dropbox (ASU)/"
droppath = "Projects/DCA/GenSeqs/"
# fullpath = macpath + droppath
fullpath = upath + droppath


pdbid = 'rnafam'
pdbfile = fullpath + "rnafam.pdb"
ATOMS = []
#EDITED TO Insert RNA P Chain
structure = Bio.PDB.PDBParser().get_structure(pdbid, pdbfile)
model = structure[0]
for residue in model.get_residues():
    tags= residue.get_full_id()
    if tags[3][0] == " ":
        atoms=residue.get_atoms()
        for atom in atoms:
            if atom.get_id() == 'CA':
                ca = atom
                ATOMS.append(ca)
for Atom in model.get_atoms():
    if Atom.get_id() == 'P':
        p = Atom
        ATOMS.append(p)

dists=[]
for atom in ATOMS:
    for ATOM in ATOMS:
        if atom.get_id()=='CA' and ATOM.get_id()=='P':
            xi, yi, zi = atom.get_coord()
            xj, yj, zj = ATOM.get_coord()
            dist = math.sqrt(((xi-xj)**2+(yi-yj)**2+(zi-zj)**2))
            dists.append((atom.get_parent().get_parent().get_id(),atom.get_parent().get_id()[1], ATOM.get_parent().get_id()[1],dist))

for c,x,y,d in dists:
    if d< 10:
        print(c,x,y,d)

# SEQ FROM PDB FILE
# CGCCTAGGTTGGGTAGGGTGGTGGCC