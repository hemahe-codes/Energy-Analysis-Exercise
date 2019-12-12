from sympy import *
from Bio.PDB import *
from optparse import OptionParser
from Bio.PDB.NACCESS import *
import os
from Bio.PDB.PDBParser import PDBParser
from residue_library import ResiduesDataLib
from forcefield import VdwParamset
from Bio.PDB.NeighborSearch import NeighborSearch
import matplotlib.pyplot as plt
import numpy as np

def solvation_1_Atom(at):
    ''' this function computes the solvation for one atom'''
    if at.element !='H':
        sigma = at.xtra['vdw'].fsrf
        surface = float(at.xtra['EXP_NACCESS'])
        return sigma * surface
    else:
        return 0


def vdw_interaction(at1,at2):
    '''this function computes the Vdw interaction for two atoms'''
    distance = at1-at2
    if distance > 0:
        try:
            epsilon = (np.sqrt((at1.xtra['vdw'].eps)*(at2.xtra['vdw'].eps)) )
            sigma_at1 = at1.xtra['vdw'].rvdw
            sigma_at2 = at2.xtra['vdw'].rvdw
            distance = at1-at2
            return  4*epsilon*(((sigma_at1/distance)**12 - (sigma_at2/distance)**6))
            
        except:
            return 0
    return 0

        


def electrostatic_interaction(at1,at2):
    ''' this function computes the electrostatic interaction for two atoms '''
    result = 0
    distance = at1-at2
    if distance > 0:
        try:
            charge_at1 = at1.xtra['charge']
            charge_at2 = at2.xtra['charge']
            mehler_solmajer_dielectric = (86.9525 / (1 - 7.7839* np.exp(-0.3153-distance))) - 8.5525
            return (332.16*((charge_at1*charge_at2)/(mehler_solmajer_dielectric*distance)))
        except:
            return 0
    return 0   

def get_contacts(at_chain_a,all_atoms,verbose,max_dist):
    final=[]
    progress = 0
    ns = NeighborSearch(all_atoms)
    for at1, at2 in ns.search_all(max_dist+1.5):
        progress+=1
        if len(verbose)>0:
            print(verbose,progress) 
        if (at1.get_full_id()[2] == 'A' and at2.get_full_id()[2] == 'B') or (at2.get_full_id()[2] == 'A' and at1.get_full_id()[2] == 'B'):
            final.append(at1)
            final.append(at2)
    return final 

			

usage = "USAGE: python GetContacts.py --f1 FirstPDB --f2 SeconPDB --c1 FirstPDBChains --c2 SecondMoleculeChains [--c ContactCutoff] [--i InterfaceCutoff] [--o1 output1] [--o2 output2] \n"
parser = OptionParser(usage=usage)

parser.add_option("--f1",help=" molecule pdb", dest="f1")
parser.add_option("--c1",help=" molecule chains", dest="c1")
parser.add_option('--res',help= 'Residue library', dest='residue_library')
parser.add_option('--vdw',help= 'Vdw', dest='vdw')
parser.add_option("--sn",help="structure name", dest="sn")

(options, args) = parser.parse_args()

print('Getting structure...')

residue_library = ResiduesDataLib(options.residue_library)
ff_params = VdwParamset(options.vdw)
str_1 = PDBParser().get_structure('molecule', options.f1) 


structure_atoms = Selection.unfold_entities(str_1,'A')

srf = NACCESS_atomic(str_1[0],naccess_binary ='/home/helena/Desktop/Biophysics_seminar_5/Biophysics/soft/NACCESS/naccess' )

max_dist = 8

chains_1 = options.c1


atoms_1 = Selection.unfold_entities(str_1, 'C') # C for chains


print('Getting atoms from chain A...')
at_chain_b = atoms_1[1].get_atoms()

print('Getting atoms from chain B...')
at_chain_a = atoms_1[0].get_atoms()
all_atoms = str_1.get_atoms()
atoms_str = []
for i in all_atoms:
    atoms_str.append(i)


print('Starting Neighbour Search...')

interface = get_contacts(at_chain_a,atoms_str,"First molecule, neighbour search of atom ",5)
interface_residues = Selection.unfold_entities(interface, 'R')

# To save interface residues from each chain in different text files:

A=[]
B=[]
for atom in interface:
    if atom.get_full_id()[2] == 'A':
        if atom.get_parent().get_id()[1] not in A:
            A.append(atom.get_parent().get_id()[1])
    if atom.get_full_id()[2] == 'B':
        if atom.get_parent().get_id()[1] not in B:
            B.append(atom.get_parent().get_id()[1])

with open('chainA_%s.txt' % (options.sn), 'w') as file:
            for residue in A:
                file.write('%s%s' % (residue,' ' ))
with open('chainB_%s.txt' % (options.sn), 'w') as file:
            for residue in B:
                file.write('%s%s' % (residue,' ' ))



for at in structure_atoms:
    resname = at.get_parent().get_resname()
    params = residue_library.get_params(resname, at.id)
    if not params:
        print("ERROR: residue/atom pair not in library (" + resname + ' ' + at.id + ')')
        continue
    at.xtra['atom_type'] = params.at_type
    at.xtra['charge'] = params.charge
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    at.xtra['solvation'] = solvation_1_Atom(at)
 
 # To compute energies

solvation_list = []
electro = []
vdw = []
energy = [1]
energy= energy*len(interface_residues)

for residue in interface_residues:
    solvation_residue =[]
    vdw_residue = []
    electro_residue = []
    for atom in residue:
        solvation_residue.append(solvation_1_Atom(atom))
        for at in interface:
            vdw_residue.append(vdw_interaction(atom,at))
            electro_residue.append(electrostatic_interaction(atom,at))
    print(vdw_residue)
    solvation_list.append(sum(solvation_residue))
    electro.append(sum(electro_residue))
    vdw.append(sum(vdw_residue))
                


for i in range(len(interface_residues)):
    energy[i] = solvation_list[i] + electro[i] + vdw[i]



# Determining the effect of each residue:

copy_energies = energy.copy()
x = {}

for i in range(len(energy)):
    copy_energies[i] = 0
    global_energy = 0
    for k in copy_energies:
        global_energy = global_energy + k
    x[i] = global_energy
    copy_energies = energy.copy()


# Plot the results:

plt.plot(list(x.keys()),list(x.values()))
plt.show()
        
