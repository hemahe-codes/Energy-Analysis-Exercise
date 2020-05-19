#! /usr/bin/python3

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import os
import matplotlib.pyplot as plt
import math
import sys

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.NeighborSearch import NeighborSearch

from residue_library import ResiduesDataLib


NACCESS_BINARY = 'path to naccess'

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

class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self, data):
        self .id = data[0] + ':' + data[1]
        self.at_type = data[2]
        self.charge = float(data[3])
class VdwParamset():
    def __init__(self, file_name):
        self.at_types = {}
        try:
            fh = open(file_name, "r")
        except OSError:
            print("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.at_types[data[0]] = AtomType(data)
        self.ntypes = len(self.at_types)
        fh.close()

class AtomType():
    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
  
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)

def add_atom_parameters(st, res_lib, ff_params):
    ''' Adds parameters from libraries to atom .xtra field
        For not recognized atoms, issues a warning and put default parameters
    '''
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = res_lib.get_params(resname, at.id)
        if not params:
            print("WARNING: residue/atom pair not in library ("+atom_id(at) + ')')
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

def get_interface(st, dist):
    ''' Detects interface residues within a distance(dist)
        Assumes two chains, i.e. a unique interface set per chain.
    '''
    select_ats = []
    for at in st.get_atoms():
        # Skip Hydrogens to reduce time
        if at.element != 'H':
            select_ats.append(at)
    nbsearch = NeighborSearch(select_ats)
    interface = {}
    # Sets are more efficient than lists. Use sets when order is not relevant
    for ch in st[0]:
        interface[ch.id] = set()

    for at1, at2 in nbsearch.search_all(dist):
        #Only different chains
        res1 = at1.get_parent()
        ch1 = res1.get_parent()
        res2 = at2.get_parent()
        ch2 = res2.get_parent()
        if ch1 != ch2:
            interface[ch1.id].add(res1)
            interface[ch2.id].add(res2)
    return interface

def calc_int_energies(st, res):
    '''Returns interaction energies (residue against other chains)
        for all atoms and for Ala atoms
    '''
    elec = 0.
    elec_ala = 0.
    vdw = 0.
    vdw_ala = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    elec_ala += e
                e = vdw_int(at1, at2, r)
                vdw += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    vdw_ala += e
    return elec, elec_ala, vdw, vdw_ala

def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    '''Electrostatic interaction energy between two atoms at r distance'''
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

def calc_solvation(st, res):
    '''Solvation energy based on ASA'''
    solv = 0.
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala       
        
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

add_atom_parameters(st, residue_library, ff_params)

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

counter = 0
print(counter)

for ch in st[0]:
    for res in ch.get_residues():
        counter=+1
        print(counter)
        if args.cutoff_dist > 0 and res not in interface[ch.id]:
            continue
        elec[res], elec_ala[res], vdw[res], vdw_ala[res] = calc_int_energies(st[0], res)
        solvAB[res], solvAB_ala[res] = calc_solvation(st[0], res)
        solvA[res], solvA_ala[res] = calc_solvation(
            st_chains[ch.id],
            st_chains[ch.id][0][ch.id][res.id[1]]
        )
        totalIntElec += elec[res]
        totalIntVdw += vdw[res]
        totalSolv += solvAB[res]
        totalSolvMon[ch.id] += solvA[res]
        total += elec[res] + vdw[res] + solvAB[res] - solvA[res]
print("Interaction energy based in interface residues only")
print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
print('{:20}: {:11.4f}'.format('DGintAB-A-B', total))
print("\#R {} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}".format(args.pdb_file.name, totalIntElec, totalIntVdw, totalSolv, totalSolvMon[chids[0]], totalSolvMon[chids[1]], total))
#Ala scanning not required. Finishing here
import sys
sys.exit()

print("Ala Scanning: DDGs for X->Ala mutations on interface residues")
for ch in st[0]:
    for res in ch.get_residues():
        if args.cutoff_dist > 0 and res not in interface[ch.id]:
            continue
        print(
            '{:11} {:11.4f}{:11.4f}{:11.4f}{:11.4f}{:11.4f}'.format(
                en.residue_id(res),
                - elec[res] + elec_ala[res],
                - vdw[res] + vdw_ala[res],
                - solvAB[res] + solvAB_ala[res],
                - solvA[res] + solvA_ala[res],
                - elec[res] + elec_ala[res] - vdw[res] + vdw_ala[res] -solvAB[res] +\
                    solvAB_ala[res] -solvA[res] + solvA_ala[res]
            )
        )
