import os
import sys
import itertools


cleave_pro = True
min_length = 7
max_mass = 4600
max_mod = 5
max_miscleavage = 2

debug_verbose = True
def debug_print(msg):
    if debug_verbose:
        print(msg)

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


i=0
def fasta(f):
    data = {}

    with open(f) as lines:
        title, seq = "", ""
        for l in lines:
            l = l.strip()
            global i
            if l.startswith(">"):
                title = l
                seq = ""

                if i > 1000000:
                    break
                i+=1
                continue

            seq += l
            data[title] = seq

    return data


proteins = fasta("Human_sp_iso_2019-08-27.fasta")
print(i, len(proteins))


fragments = set()
fragments_acetyl = set()

# TODO what do we do for "X" ?
allseq = proteins.values()
allseq = set(filter(lambda x: "X" not in x, allseq))
del proteins

# DEBUG allseq = ["HSFJKHAJKFSHJKFSAHJKCHSJKCAJSJRWKJEKDLSJKRNWKENKLWQNDQLKFNK"]
for s in allseq:
    #print(s)

    # Cleave at Lysine(Lys, K) and Arginine (Arg, R)
    csite = [-1]
    for i, aa in enumerate(s):
        if aa in "KR":
            # TODO Not taking account in proline
            csite.append(i)
    if len(s)-1 not in csite:
        csite.append(len(s)-1)

    for i, cs1 in enumerate(csite):
        for j, cs2 in enumerate(csite):
            if j <= i:
                continue
            if j-i-1 > max_miscleavage:
                continue

            frag = s[cs1+1:cs2+1]
            if len(frag) < min_length:
                continue
            #print(" "*(cs1+1) + frag)
            fragments.add(frag)

            if cs1 == -1:
                fragments_acetyl.add(frag)




print(len(allseq))
print(len(fragments))

del allseq


outputlines = []
for f in fragments:
    possible_mods = []
    # debug_print(f)

    base_mass = 2*atom_mass["H"] + atom_mass["O"]
    for i, aa in enumerate(f):
        if aa == "M":
            possible_mods.append((i+1, "Oxidation"))
        base_mass += aa_mass[aa]

    if base_mass > max_mass:
        continue

    if f in fragments_acetyl:
        possible_mods.append((0, "Acetyl"))

    for numOfMod in range(min(len(possible_mods)+1, max_mod+1)):
        selected_mods = itertools.combinations(possible_mods, numOfMod)
        for selected_mod in selected_mods:
            #debug_print(selected_mod)

            selected_mod = list(selected_mod)
            modstr_almost = []

            for i, aa in enumerate(f):
                if aa == "C":
                    selected_mod.append((i+1, "Carbamidomethyl"))

                if aa == "K":
                    selected_mod.append((i+1, "TMT"))

            if (0, "Acetyl") not in selected_mod:
                selected_mod.append((0, "TMT"))

            newmass = base_mass
            selected_mod.sort()
            for tmp in selected_mod:
                newmass += mod_mass[tmp[1]]
                modstr_almost.append("%s|%s" % tmp)

            if newmass > max_mass:
                continue

            l = [f, "|".join(modstr_almost), "%0.2f" % newmass]
            outputlines.append(",".join(l) + "\n")





outputfile = open("../../delivery/human_library.csv", "w")
outputfile.write("seq,modifications\n")
for l in sorted(outputlines):
    outputfile.write(l)







