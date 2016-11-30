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
import subprocess
#from prep_calc import *
from geometry import *
from atom3D import *
from globalvars import globalvars
from mol3D import*



########### UNIT CONVERSION
HF_to_Kcal_mol = 627.503
#def maximum_ML_dist(mol):
#    core = mol.getAtom(mol.findMetal()).coords()
#    max_dist = 0
#    for atom_inds in mol.getBondedAtoms(mol.findMetal()):
#        dist = distance(core,mol.getAtom(atom_inds).coords())
#        if (dist > max_dist):
#            max_dist = dist
#    return max_dist

#def minimum_ML_dist(mol):
#    core = mol.getAtom(mol.findMetal()).coords()
#    min_dist = 1000
#    for atom_inds in mol.getBondedAtoms(mol.findMetal()):
#        dist = distance(core,mol.getAtom(atom_inds).coords())
#        if (dist < min_dist) and (dist > 0):
#            min_dist = dist
#    return min_dist

def translate_job_name(job):
    base = os.path.basename(job)
    base = base.strip("\n")
    base_name = base.strip(".in")
    base_name = base_name.strip(".done")
    print(base_name)
    low_name = str(base_name).lower()
    ll = (str(base_name)).split("_")
    ID = int(ll[0])
    metal = ll[1].lower()
    ox = int(ll[2])
    eqlig = int(ll[4])

    axlig1 = int(ll[6])
    axlig2 = int(ll[8])
    spin = int(ll[9])
    metal_spin_dictionary = {'co':{2:[2,4],3:[1,5]},
                              'cr':{2:[3,5],3:[2,4]},
                              'fe':{2:[1,5],3:[2,6]},
                              'mn':{2:[2,6],3:[3,5]},
                              'ni':{2:[1,3]}}
    these_states = metal_spin_dictionary[metal][ox]
    if spin == these_states[0]:
            spin_cat = 'LS'
    elif spin == these_states[1]:
            spin_cat = 'HS'
    else:
        print('critical erorr, unknown spin: '+ str(spin))
    return ID,low_name,base_name,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat


class DFTRun:
    """ This is a class for each run"""
    numRuns = 0
    def __init__(self,name):
        self.numRuns += 1
        self.name = name
        self.outpath  = 'undef'
        self.geopath = 'undef'
        self.geo_exists = False
        self.output_exists = False
        self.converged = False
        self.time = 'undef'
        self.energy = 0
        self.spin = 'undef'
        self.bound  = False
        self.lig1 = 'undef'
        self.lig2 = 'undef'
        self.rmsd = 'undef'
        self.spin_cat = 'undef'
        self.ss_act = 0
        self.ss_target = 0
    def extract_geo(self):
            cmd_str = 'python optgeo_extract.py '+ self.scrpath + ' ' + self.optgeopath
            print(cmd_str)
            subprocess.call(cmd_str,shell=True)
    def extract_prog(self):
            cmd_str = 'python optgeo_extract.py '+ self.scrpath + ' ' + self.progpath
            print(cmd_str)
            subprocess.call(cmd_str,shell=True)

    def obtain_mol3d(self):
        this_mol = mol3D()
        if os.path.exists(self.geopath):
                this_mol.readfromxyz(self.geopath)
        elif  os.path.exists(self.progpath):
                this_mol.readfromxyz(self.progpath)
        self.mol = this_mol
        init_mol = mol3D()
        init_mol.readfromxyz(self.init_geopath)
        self.init_mol = init_mol
    def obtain_rmsd(self):
        self.rmsd = self.mol.rmsd(self.init_mol)
        print('rmsd = '+ str(self.rmsd))
#    def obtain_ML_dists(self):
#        self.min_dist = minimum_ML_dist(self.mol)
#        self.max_dist = maximum_ML_dist(self.mol)
    def configure(self,metal,ox,spin,has_o2,eqlig,axlig1,axlig2):
        self.metal = metal
        self.ox = ox
        self.spin = spin
        self.bound = has_o2
        self.lig1 ='cl'
        self.lig2 = 'cl'

class Comp:
    """ This is a class for each unique composition and configuration"""
    def __init__(self,name):
        self.name = name
        self.ox =' undef'
        self.metal= 'undef'
        self.lig1 = 'undef'
        self.lig2 = 'undef'
        self.LSenergy = 'undef'
        self.HSenergy= 'undef' 
        self.HSss_target = 'undef'
        self.LSss_target = 'undef'
        self.HSss_act = 'undef'
        self.LSss_act = 'undef'
        self.LS_rmsd = 'undef'
        self.HS_rmsd = 'undef'


def writeprops(extrct_props,newfile,do_strip):
    string_to_write = ','.join([str(word) for word in extrct_props ])
    newfile.write(string_to_write)
    newfile.write("\n")
    return 
def scfextract(a_run,list_of_props):
    extrct_props = []
    for keys in list_of_props:
        extrct_props.append(a_run.__dict__[str(keys)])
    return extrct_props
# pass the name of the species/test to the script as only agrument
def test_terachem_go_convergence(job):
    ### get paths
    basic_path  = '/home/jp/solpm7test/d0/'
    ID,low_name,base_name,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat = translate_job_name(job)
    ### flag
    converged =  False
    geo_exists = False
    ### test if geo exits
    this_run=DFTRun(low_name)
    this_run.geopath = (basic_path + 'optimized_geo/' + low_name + ".xyz")
    this_run.progpath = (basic_path + 'prog_geo/' + low_name + ".xyz")
    this_run.init_geopath = (basic_path + 'initial_geo/' + base_name + ".xyz")
    this_run.scrpath = (basic_path + 'scr/' + low_name +"/optim.xyz")
    this_run.outpath = (basic_path + 'outfiles/' + base_name + ".out")
    this_run.configure(metal,ox,spin,False,eqlig,axlig1,axlig2)
    this_run.ID = ID
    this_run.converged = False
    this_run.spin_cat = spin_cat
    print('run is set up')
    if os.path.exists(this_run.geopath):
        this_run.geo_exists = True
 #       this_run.obtain_mol3d()
    if os.path.exists(this_run.outpath):
        ### file is found,d check if converged
        with open(this_run.outpath) as f: 
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('Converged!') != -1:
                    found_conv = True
    #                print(lines)
                if str(lines).find('Optimization Converged.') != -1:
                   found_conv = True
                   print(lines)
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy =str(lines.split()[2])
                    found_data = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time=str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str=(lines.split())
                    this_run.ss_act =float( this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
            print('found all')
        if this_run.converged:
                print('run converged')
                if this_run.geo_exists:
                        print('geo exists for ' + this_run.name)
                if not this_run.geo_exists:
                        try:
                                this_run.extract_geo()
                        except:
                                print("ERROR: scr not found for" + str(this_run.geopath))
                if os.path.exists(this_run.geopath):
                        if os.path.exists(this_run.init_geopath):
                                print('both paths exist')
                                this_run.obtain_mol3d()
                                try:
                                        this_run.obtain_rmsd()
                                except:
                                        this_run.rmsd = "undef"
        if not this_run.converged:
                print(' \n job  ' + str(this_run.outpath) + ' not converged\n')
                try:
                        this_run.extract_prog()
                        this_run.obtain_mol3d()
                        try:
                            this_run.obtain_rmsd()
                        except:
                            this_run.rmsd = "undef"

                except:
                        print("ERROR: scr not found for" + str(this_run.progpath))
    print('run is returned')
    return this_run
def test_terachem_sp_convergence(job):
    ### get paths
    path_dictionary = setup_paths()
    gen,slot,gene,spin,base_name = translate_job_name(job)
    ### flag
    converged =  False
    ### test if geo exits
    this_run=DFTRun(base_name)
    this_run.outpath = (path_dictionary["out_path" ] + "/gen_" + str(gen) +"/"
                           + base_name + ".out")
    print("checking ",this_run.outpath)
    this_run.spin = spin
    this_run.gene =  gene
    if os.path.exists(this_run.outpath):
        ### file is found,d check if converged
        with open(this_run.outpath) as f: 
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('Running Mulliken') != -1:
                    found_conv = True
                if str(lines).find('FINAL ENERGY') != -1:
                    this_run.energy =str(lines.split()[2])
                    found_data = True
                if str(lines).find('Total processing time') != -1:
                    this_run.time=str(lines.split()[3])
                    found_time = True
                if str(lines).find('SPIN S-SQUARED') != -1:
                    this_str=(lines.split())
                    this_run.ssq =float( this_str[2])
                    this_run.star = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
            print('run converged : ' + str(this_run.name))
    return this_run


def process_runs(LS_runs,HS_runs):
    final_results=dict()
    matched = False
    number_of_matches  = 0
    print('processing all converged runs')
    for ID in LS_runs.keys():
        matched = 0 
        LS_run = LS_runs[ID]
        this_name = LS_run.name
        this_ID = ID
        if ID in HS_runs.keys():
            HS_run = HS_runs[ID]
            matched = True
            number_of_matches += 1
        if matched:
            print('matched ID: '+ str( ID) + ' files ' + str(HS_run.name) + ' and ' + str(LS_run.name))
            final_results[this_ID] = Comp(this_ID)
            final_results[this_ID].ID = this_ID
            final_results[this_ID].LSenergy = str(float(LS_run.energy)*HF_to_Kcal_mol)
            final_results[this_ID].HSenergy = str(float(HS_run.energy)*HF_to_Kcal_mol)
            final_results[this_ID].splitenergy = str((float(HS_run.energy)*HF_to_Kcal_mol)-(float(LS_run.energy)*HF_to_Kcal_mol))
#            final_results[this_ID].ligno1 = LS_run.ligno1
#            final_results[this_ID].ligno2 = LS_run.ligno2
            final_results[this_ID].lig1 = LS_run.lig1
            final_results[this_ID].lig2 = LS_run.lig2
            final_results[this_ID].HSss_act = HS_run.ss_act
            final_results[this_ID].LSss_act = LS_run.ss_act
            final_results[this_ID].LSss_target = LS_run.ss_target
            final_results[this_ID].HSss_target = HS_run.ss_target
            final_results[this_ID].LS_rmsd = LS_run.rmsd
            final_results[this_ID].HS_rmsd = HS_run.rmsd
 
        else:
            print('unmatched ID: '+ str( ID) + ' files ' + str(LS_run.name)+ ' has no partner' )

    return final_results

full_list_of_props = ['ID','name','converged','spin','lig1','lig2','ss_act','ss_target','energy','rmsd']
summary_list_of_props = ['ID','LSenergy','HSenergy','LSss_act','LSss_target','HSss_act','HSss_target','LS_rmsd','HS_rmsd','splitenergy']
new_file = open('full_results.txt','w')
new_file.write(','.join(full_list_of_props) + '\n')
summary_file = open('summary_results.txt','w')
summary_file.write(','.join(summary_list_of_props) + '\n')
targetpaths=sorted(glob.glob("completejobs/*.in"))
LS_runs = dict()
HS_runs = dict()
for files in targetpaths:
        print(files)
        this_run = test_terachem_go_convergence(files)
        extrct_props = scfextract(this_run,full_list_of_props)
        writeprops(extrct_props,new_file,1)
        if this_run.converged:
           if this_run.spin_cat == 'LS':
                LS_runs.update({this_run.ID:this_run})
           else:
                HS_runs.update({this_run.ID:this_run})
targetpaths=sorted(glob.glob("jobs/*.in"))
print('****************************************************************')
print('****************************************************************')

for files in targetpaths:
        print('target  is : '+ files)
        this_run = test_terachem_go_convergence(files)
        print('run returned to loop')
        extrct_props = scfextract(this_run,full_list_of_props)
        print('writing: '+ str(extrct_props))
        writeprops(extrct_props,new_file,1)
        print('written')
        if this_run.converged:
           if this_run.spin_cat == 'LS':
                LS_runs.update({this_run.ID:this_run})
           else:
                HS_runs.update({this_run.ID:this_run})

print('****************************************************************')
print('****************************************************************')

final_results = process_runs(LS_runs,HS_runs)
for keys in sorted(final_results):
    this_result = final_results[keys]
    extrct_props = scfextract(this_result,summary_list_of_props)
    writeprops(extrct_props,summary_file,1)



new_file.close()
summary_file.close()
