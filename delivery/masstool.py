import os
import sys

# Residue mass...
aa_mass = {}
aa_table = """glycine	G	GLY	C2H3NO	57.021463735	57.05132
alanine	A	ALA	C3H5NO	71.037113805	71.0779
serine	S	SER	C3H5NO2	87.032028435	87.0773
proline	P	PRO	C5H7NO	97.052763875	97.11518
valine	V	VAL	C5H9NO	99.068413945	99.13106
threonine	T	THR	C4H7NO2	101.047678505	101.10388
cysteine	C	CYS	C3H5NOS	103.009184505	103.1429
leucine	L	LEU	C6H11NO	113.084064015	113.15764
isoleucine	I	ILE	C6H11NO	113.084064015	113.15764
asparagine	N	ASN	C4H6N2O2	114.042927470	114.10264
aspartic_acid	D	ASP	C4H5NO3	115.026943065	115.0874
glutamine	Q	GLN	C5H8N2O2	128.058577540	128.12922
lysine	K	LYS	C6H12N2O	128.094963050	128.17228
glutamic_acid	E	GLU	C5H7NO3	129.042593135	129.11398
methionine	M	MET	C5H9NOS	131.040484645	131.19606
histidine	H	HIS	C6H7N3O	137.058911875	137.13928
phenylalanine	F	PHE	C9H9NO	147.068413945	147.17386
selenocysteine	U	SEC	C3H5NOSe	150.953633405	150.3079
arginine	R	ARG	C6H12N4O	156.101111050	156.18568
tyrosine	Y	TYR	C9H9NO2	163.063328575	163.17326
tryptophan	W	TRP	C11H10N2O	186.079312980	186.2099
pyrrolysine	O	PYL	C12H19N3O2	237.147726925	237.29816"""
for l in aa_table.strip().split("\n"):
    lname, l1, l3, comp, mass1, mass2 = l.split()
    aa_mass[l1] = float(mass1)

atom_mass = {}
atom_str = """Hydrogen	H	1.007825035	1.00794
Carbon	C	12.0000000	12.0107
Nitrogen	N	14.003074	14.0067
Oxygen	O	15.99491463	15.9994
Phosphorus	P	30.973762	30.973761
Sulphur	S	31.9720707	32.065
Selenium	Se	79.9165196	79.96
proton	proton	1.00727646688	1.00727646688"""
for l in atom_str.strip().split("\n"):
    lname, abr, mass1, mass2 = l.split()
    atom_mass[abr] = float(mass1)

mod_mass = {
    "TMT": 12*atom_mass["C"] + 20*atom_mass["H"] + 2*atom_mass["O"] + 2*atom_mass["N"],
    "Oxidation": atom_mass["O"],
    "Carbamidomethyl": 2*atom_mass["C"] + 3*atom_mass["H"] + 1*atom_mass["O"] + 1*atom_mass["N"],
    "Acetyl": 2*atom_mass["C"] + 2*atom_mass["H"] + 1*atom_mass["O"],
}




def peptide_mw(pep):
    res = atom_mass["H"]*2 + atom_mass["O"]
    for aa in pep:
        res += aa_mass[aa]
    return res


if __name__ == "__main__":
    pep = "PLHKYPVWLWK"
    dmw = 2138.28256686736
    mw = peptide_mw(pep) #+ mod_mass["Carbamidomethyl"] #+ mod_mass["Acetyl"] + mod_mass["Oxidation"]
    print(mw)
    print(dmw - mw)

    print(mod_mass["TMT"] )






