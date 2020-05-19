#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import os
import matplotlib.pyplot as plt

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select

from residue_library import ResiduesDataLib
from forcefield import VdwParamset
import energies as en

NACCESS_BINARY = '/home/helena/Escritorio/Biophysics-master_2/Biophysics-master/soft/NACCESS/naccess'

parse_cmd = argparse.ArgumentParser(
    prog='binding',
    description='binding energy calculation'
)

parse_cmd.add_argument(
    '--rlib',
    action='store',
    dest='reslib_file',
    default='data/aaLib.lib',
    help='Residue Library'
)
parse_cmd.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default='data/vdwprm',
    help='Vdw parameters'
)

parse_cmd.add_argument(
    '--dist',
    action='store',
    dest='cutoff_dist',
    default=8.0,
    type=float,
    help='Cutoff distance for determining the interface (0: use all residues):'
)
parse_cmd.add_argument('pdb_file', help='Input PDB', type=open)

args = parse_cmd.parse_args()

print("PDB.filename:", args.pdb_file.name)
print("Residue Lib.:", args.reslib_file)
print("PDB.filename:", args.vdwprm_file)
print("Distance:", args.cutoff_dist)

# Loading Libraries
# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib(args.reslib_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing', args.pdb_file)
# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# assign data types, and charges from libraries
# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Possible errors on N-term and C-Term atoms
# Possible errors on HIS alternative forms

en.add_atom_parameters(st, residue_library, ff_params)

# Calculating surfaces
# The specific PATH to naccess script (in soft) is needed
# ASA goes to .xtra field directly

srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

# Prepare surfaces for the separate chains
# Alternatively the twp PDB files can be prepared outside and parsed here

io = PDBIO()
st_chains = {}
# Using BioIO trick (see tutorial) to select chains
class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid

    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0

for ch in st[0]:
    io.set_structure(st)
    io.save('tmp.pdb', SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
    en.add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')

## Interface residues
if args.cutoff_dist > 0.:
    interface = en.get_interface(st, args.cutoff_dist)

## Initiatlize Energy aggregates
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

totalIntElec = 0.
totalIntVdw = 0.
totalSolv = 0.
totalSolvMon = {}
## We get the chsin ids,not always they are A and B
chids = []
for ch in st[0]:
    chids.append(ch.id)
    totalSolvMon[ch.id] = 0

total = 0.

for ch in st[0]:
    for res in ch.get_residues():
        if args.cutoff_dist > 0 and res not in interface[ch.id]:
            continue
        elec[res], elec_ala[res], vdw[res], vdw_ala[res] = en.calc_int_energies(st[0], res)
        solvAB[res], solvAB_ala[res] = en.calc_solvation(st[0], res)
        try:
            solvA[res], solvA_ala[res] = en.calc_solvation(st_chains[ch.id],st_chains[ch.id][0][ch.id][res.id[1]])
        except:
            solvA[res],solvA_ala[res] = 0,0
        totalIntElec += elec[res]
        totalIntVdw += vdw[res]
        totalSolv += solvAB[res]
        totalSolvMon[ch.id] += solvA[res]
        total += elec[res] + vdw[res] + solvAB[res] - solvA[res]
#Ala scanning not required. Finishing here
import sys
#sys.exit()

with open('chainA_%s.txt' % ('4NCO'), 'w') as file:
    x=[]
    y=[]
    for ch in st[0]:
        for res in ch.get_residues():
            if args.cutoff_dist > 0 and res not in interface[ch.id]:
                continue
            if (elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +solvAB_ala[res] -solvA[res] + solvA_ala[res]) > 8.214602:
                file.write('%s%s' % (en.residue_id(res)[5:],' ' ))
            x.append(en.residue_id(res))
            y.append(elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +solvAB_ala[res] -solvA[res] + solvA_ala[res])
    print(max(y))
    print(min(y))
    plt.plot(x,y)
    plt.show()
