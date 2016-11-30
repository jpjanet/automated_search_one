import glob
import string
import sys
import os
import numpy as np
import math
import random
import string
import numpy
import pybel
from geometry import *
from atom3D import *
from globalvars import globalvars
from mol3D import*
#name    metal   ox  axlig_charge    eqlig charge    axlig_dent  eqlig_dent  axlig_connect   eqlig_connect   axlig_natoms    eqlig_natoms    axlig_mdelen    eqlig_mdelen


########### UNIT CONVERSION
HF_to_Kcal_mol = 627.503
def maximum_ML_dist(mol):
    core = mol.getAtom(mol.findMetal()).coords()
    max_dist = 0
    for atom_inds in mol.getBondedAtoms(mol.findMetal()):
        dist = distance(core,mol.getAtom(atom_inds).coords())
        if (dist > max_dist):
            max_dist = dist
    return max_dist

def minimum_ML_dist(mol):
    core = mol.getAtom(mol.findMetal()).coords()
    min_dist = 1000
    for atom_inds in mol.getBondedAtoms(mol.findMetal()):
        dist = distance(core,mol.getAtom(atom_inds).coords())
        if (dist < min_dist) and (dist > 0):
            min_dist = dist
    return min_dist

def compare_tan_to_ligands(path):
    gl = glob.glob("Ligands/*.xyz")
    this_lig = mol3D()
    ref_OBmol = this_lig.getOBmol(path,"xyzf")
    best_tan = 0
    best_lig = "nothing !"
    for ligand_paths  in gl:
        this_lig = mol3D()
        this_OBmol = this_lig.getOBmol(ligand_paths,"xyzf")
        this_tan = ref_OBmol.calcfp() | this_OBmol.calcfp()
        if this_tan > best_tan:
            best_tan = this_tan
            best_lig = ligand_paths
            print("ligand is " + ligand_paths)
            print("tanimoto inedx : " + str(this_tan) )
    return best_tan,best_lig
class octahedral_complex:
    def __init__(self,name):
       self.path = "undef"
       self.mol = mol3D() 
       self.liglist = list()
       self.ligdents = list()
       self.ligcons = list()
    def obtain_mol3d(self):
        this_mol = mol3D()
        this_mol.readfromxyz(self.path)
        self.mol = this_mol
    def ligand_breakdown(self):
        metal_index = self.mol.findMetal()
#        print(metal_index)
        bondedatoms = self.mol.getBondedAtoms(metal_index)
        if ((self.path =="optimized_geo/test_29_20.xyz") or (self.path == "optimized_geo/test_29_00.xyz") or
            (self.path =="optimized_geo/test_30_20.xyz") or (self.path == "optimized_geo/test_30_00.xyz")):
            correct = list()
            for ba in bondedatoms:
                if (self.mol.getAtom(ba).symbol() == "O"):
                    correct.append(ba)
            bondedatoms = correct
        if ((self.path =="optimized_geo/test_39_20.xyz") or (self.path == "optimized_geo/test_39_00.xyz")):
#            (self.path =="optimized_geo/test_30_20.xyz") or (self.path == "optimized_geo/test_30_00.xyz")):
#            list_of_o = [ atom for atom in range(0,self.mol.natoms()) if (self.mol.getAtom(atom).symbol() == "O")]
            list_of_o = [ atom for atom in range(0,self.mol.natoms) if (self.mol.getAtom(atom).symbol() == "O")]

            bondedatoms += list_of_o

        if  (self.path == "optimized_geo/test_27_00.xyz"):
            list_of_o = [ atom for atom in range(0,self.mol.natoms) if (self.mol.getAtom(atom).symbol() == "O")]

            bondedatoms = list_of_o


        counter = 0
        for atom in bondedatoms:
                print('this atom type is ' + self.mol.getAtom(atom).symbol())
                print('conection number ' + str(atom) + " of " + str(bondedatoms))
                fragment = self.mol.findsubMol(atom,metal_index)
                this_cons = [x for x in fragment if (x in bondedatoms)]
                unique =  True
                for i,unique_ligands in enumerate(self.liglist):
                    print(i)
                    if fragment == unique_ligands:
                        print(fragment)
                        print("is duplicated")
#                        print("\n\n")
                        unique = False
                        matched = i
                if unique:
                   self.liglist.append(fragment)
                   self.ligdents.append(1)
                   self.ligcons.append(this_cons)
                else:
                    self.ligdents[matched] += 1
        if  (self.path == "optimized_geo/test_74_20.xyz"):
            print(bondedatoms)
            sjkasjld

    def ligand_assign(self):
        metal_index = self.mol.findMetal()
        built_ligand_list  = list()
        lig_natoms_list = list()
        n_ligs = len(self.liglist)
        max_dent = max(self.ligdents)
        min_dent = min(self.ligdents)
        print("n_ligs = " + str(n_ligs))
        print("max d = " + str(max_dent))
        print(" min_dent = " +  str(min_dent))
        for i,ligand_indices in enumerate(self.liglist):
         this_ligand = ligand(self.mol,ligand_indices,self.ligdents[i])
         this_ligand.obtain_mol3d()
         this_ligand.obtain_truncation(self.ligcons[i],3)
         built_ligand_list.append(this_ligand)

         lig_natoms_list.append(this_ligand.mol.natoms)

         if (n_ligs == 3) or (n_ligs == 4): # most common case, 
                                           # one/two equitorial and 2 axial mono
                                           # or three bidentate 

            if self.ligdents[i] == 1:  ## anything with equitorial monos will 
                                 ## have higher than 4 n_ligs
                ax_lig = i
                ax_con = self.ligcons[i]
            if (self.ligdents[i] >= 2) and (min_dent == 1):
                eq_lig = i
                eq_con = self.ligcons[i]
            if (n_ligs == 3) and (min_dent == max_dent):
                # take any 2, they are all the same
                ax_lig = 0
                eq_lig = 1
                ax_con = self.ligcons[0]
                eq_con = self.ligcons[1]
        if (n_ligs == 6): # all mono  case, 
            minz = -500
            maxz = 500
 #           print('six')
            for j,built_ligs in enumerate(built_ligand_list):
                  this_z = built_ligs.mol.centermass()[2]
            if this_z > minz:
                minz = this_z
                ax_lig = j
                ax_con = self.ligcons[j]
            if this_z < maxz:
                maxz = this_z
                not_eq = j
                allowed = range(0,6)
                allowed = [x for x in allowed  if ((x != not_eq) and (x != ax_lig))]
                eq_lig = allowed[0]
            unique_ligands = []
            ligand_counts  = list()
            ligand_records = list()

            for j,built_ligs in enumerate(built_ligand_list):
 #               print('\nHello')
                sl =  [ atom.symbol() for atom in built_ligs.mol.getAtoms()]
                unique = 1
                for i,other_sl in enumerate(unique_ligands):
                   if sorted(sl) == sorted(other_sl):
                        #duplicate
                        unique = 0
                        ligand_counts[i] +=1
                        print('\n duplicate')
                if unique == 1:
                   unique_ligands.append(sl)
                   ligand_counts.append(1)
                   ligand_records.append(j)
#            print(unique_ligands)
#            print(ligand_counts)
#            print(ligand_records)
            if (max(ligand_counts) != 4) or (min(ligand_counts) != 2):
                if (max(ligand_counts) == 6):
                   ax_lig=ligand_records[ligand_counts.index(6)]
                   eq_lig=ligand_records[ligand_counts.index(6)]
                else:
                    print('critical error, monodentates not the same')
                    print(ligand_counts)
                    print(unique_ligands)
            else:
                ax_lig=ligand_records[ligand_counts.index(2)]
                eq_lig=ligand_records[ligand_counts.index(4)]

        self.ax_ligand = built_ligand_list[ax_lig]
        self.eq_ligand = built_ligand_list[eq_lig]
        self.ax_kier = kier(self.ax_ligand.mol)
        self.ax_kier_t = kier(self.ax_ligand.trunc_mol) 
        self.eq_kier = kier(self.eq_ligand.mol)
        self.eq_kier_t = kier(self.eq_ligand.trunc_mol) 
        self.ax_natoms = lig_natoms_list[ax_lig]
        self.eq_natoms = lig_natoms_list[eq_lig]


class ligand:
    def __init__(self,master_mol,index_list,dent):
        self.master_mol  = master_mol
        self.index_list = index_list
        self.dent = dent
    def obtain_mol3d(self):
        this_mol = mol3D()
        for i in range(0,self.master_mol.natoms): 
            if i in self.index_list:
                this_mol.addAtom(self.master_mol.getAtom(i))
        self.mol = this_mol 
    def obtain_truncation(self,con_atoms,hops):
        self.trunc_mol = mol3D()
        added_list = list()
        for connections in con_atoms:
            hopped = 0
            active_set  = [connections]
            while hopped < hops:
                hopped += 1
                new_active_set = list()
                for this_atom in active_set:
                    this_atoms_neighbors =  self.master_mol.getBondedAtoms(this_atom)
                    for bound_atoms in this_atoms_neighbors:
                        if (bound_atoms in self.index_list) and (bound_atoms not in added_list):
                            self.trunc_mol.addAtom(self.master_mol.getAtom(bound_atoms))
                            added_list.append(bound_atoms)
                        [new_active_set.append(element) for element in this_atoms_neighbors]
                active_set = new_active_set

def truncate_loose(mol,nodes,hops):
        trunc_mol = mol3D()
        added_list = list()
        for node in nodes:
            hopped = 0
            active_set  = [node]
            while hopped < hops:
                hopped += 1
                new_active_set = list()
                for this_atom in active_set:
                    this_atoms_neighbors =  mol.getBondedAtoms(this_atom)
                    for bound_atoms in this_atoms_neighbors:
                        if bound_atoms not in added_list:
                            trunc_mol.addAtom(mol.getAtom(bound_atoms))
                            added_list.append(bound_atoms)
                        [new_active_set.append(element) for element in this_atoms_neighbors]
                active_set = new_active_set
        return trunc_mol



def create_graph(mol):
    index_set = range(0,mol.natoms)
    A  = numpy.matrix(numpy.zeros((mol.natoms,mol.natoms)))
    for i in index_set:
 ##       print("bonding from " + mol.getAtom(i).symbol() )
        this_bonded_atoms = mol.getBondedAtoms(i)
        for j in index_set:
            if j in this_bonded_atoms:
##                print(mol.getAtom(j).symbol() + " in bonded list")
                A[i,j] = 1
##        print("\n\n")
    return A
def remove_diagonals(matrix):
    n = matrix.shape[0]
    for i in range(0,n):
        matrix[i,i] = 0
    return matrix
def kier(mol):
    copy_mol = mol3D()
    copy_mol.copymol3D(mol)
    copy_mol.deleteHs()
    A = create_graph(copy_mol)
    n = A.shape[0]
    twopath = A*A
    remove_diagonals(twopath)
    p2 = twopath.sum()/2
    print('P2 = '+ str(p2))

    if (p2 != 0):
        two_kappa = ((numpy.power(n,3) - 5*numpy.power(n,2) + 8*n -4)
                   /(numpy.power(p2,2)))    
    else:
        two_kappa = 0 
    return(two_kappa)
def test_extract(number,alpha):
    this_complex = octahedral_complex(number)
    path = "optimized_geo/test_"+str(number) + "_" + str(alpha) +  ".xyz"
    print("\n  this is test number " + str(number))
    print(path)
    this_complex.path = path
    this_complex.obtain_mol3d()
    this_complex.ligand_breakdown()
    this_complex.ligand_assign()
    this_complex.eq_ligand.mol.writexyz("csd_ligands/test_" +str(number)+ "_eq_lig.xyz")
    this_eq_tan,this_eq_tan_match = compare_tan_to_ligands("csd_ligands/test_" +str(number)+ "_eq_lig.xyz")
    this_complex.eq_tan =  this_eq_tan
    this_complex.eq_tan_match = this_eq_tan_match
    this_complex.ax_ligand.mol.writexyz("csd_ligands/test_" +str(number)+ "_ax_lig.xyz")
    this_ax_tan,this_ax_tan_match = compare_tan_to_ligands("csd_ligands/test_" +str(number)+ "_ax_lig.xyz")
    this_complex.ax_tan =  this_ax_tan
    this_complex.ax_tan_match = this_ax_tan_match

    return this_complex
def find_kiers(name,dent):
    mol = mol3D()
    mol.readfromxyz("Ligands/" + name +".xyz")
    kappa = kier(mol)
    mol_trunc = truncate_loose(mol,range(0,dent),3)
    mol_trunc.writexyz(name + "_trunc.xyz")
    trunc_kappa = kier(mol_trunc)
    print("kappa is " + str(kappa))
    print("kappa 3 is " + str(trunc_kappa))




acac = mol3D()
acac_kappa = kier(acac)

water = mol3D()
water.readfromxyz("Ligands/water.xyz")
water_kappa = kier(water)

por = mol3D()
por.readfromxyz("Ligands/porphyrin.xyz")
por_kappa = kier(por)
por_trunc = truncate_loose(por,[0, 1, 2, 3],3)
por_trunc.writexyz("por_trunc.xyz")
por_trunc_kappa = kier(por_trunc)




pisc= mol3D()
pisc.readfromxyz("Ligands/tbuphisocy.xyz")
pisc_kappa = kier(pisc)
pisc_trunc = truncate_loose(pisc,[0],3)
pisc_trunc.writexyz("pisc_trunc.xyz")
pisc_trunc_kappa = kier(pisc_trunc)



pisc= mol3D()
pisc.readfromxyz("Ligands/tbuphisocy.xyz")
pisc_kappa = kier(pisc)
pisc_trunc = truncate_loose(pisc,[0],3)
pisc_trunc.writexyz("pisc_trunc.xyz")
pisc_trunc_kappa = kier(pisc_trunc)

sixtyfour = octahedral_complex("testing")
sixtyfour.path="optimized_geo/test_5_20.xyz"
sixtyfour.obtain_mol3d()
sixtyfour.ligand_breakdown()
sixtyfour.ligand_assign()
sixtyfour.ax_ligand.trunc_mol.writexyz("axial_trun.xyz")
sixtyfour.eq_ligand.trunc_mol.writexyz("equit_trun.xyz")
sixtyfour.eq_ligand.mol.writexyz("equit.xyz")

print("\n pisc")
find_kiers("tbuphisocy",1)
print("\n misc")
find_kiers("misc",1)
print("\n porhyrin")
find_kiers("porphyrin",4)
print("\n thiocynante")
find_kiers("thiocyanate",1)
print("\n isothiocynante")
find_kiers("isothiocyanate",1)
print("\n bipy")
find_kiers("bipy",2)
print("\n en")
find_kiers("en",2)
print("\n phen")
find_kiers("phen",2)
print("\n tbucat")
find_kiers("tbucat",2)
print("\n oxalate")
find_kiers("oxalate",2)
print("\n acac")
find_kiers("acac",2)
print("\n water")
find_kiers("water",2)


