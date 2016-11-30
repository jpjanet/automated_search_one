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
from prep_calc import *
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
        self.ss_act = 0
        self.ss_target = 0
    def extract_geo(self):
            cmd_str = 'python optgeo_extract.py '+ self.scrpath + ' ' + self.geopath
            print(cmd_str)
            subprocess.call(cmd_str,shell=True)
    def obtain_mol3d(self):
        this_mol = mol3D()
        this_mol.readfromxyz(self.geopath)
        self.mol = this_mol
    def configure(self,metal,ox,spin,has_o2,lig1,lig2):
        ligands_list  = get_ligands()
        self.metal = metal
        self.ox = ox
        self.spin = spin
        self.bound = has_o2
        self.lig1 = ligands_list[int(lig1)][0]
        self.lig2 = ligands_list[int(lig2)][0]

class Comp:
    """ This is a class for each unique composition and configuration"""
    def __init__(self,name):
        self.name = name
        self.ox =' undef'
        self.metal= 'undef'
        self.lig1 = 'undef'
        self.lig2 = 'undef'
        self.unboundenergy = 'undef'
        self.boundenergy= 'undef' 
        self.boundss_target = 'undef'
        self.unboundss_target = 'undef'
        self.boundss_act = 'undef'
        self.unboundss_act = 'undef'


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
    basic_path  = '/home/jp/minimal_models/'
    ID,low_name,base_name,metal,ox,lig1,lig2,spin,has_o2 = translate_job_name(job)
    ### flag
    converged =  False
    geo_exists = False
    ### test if geo exits
    this_run=DFTRun(low_name)
    this_run.geopath = (basic_path + 'optimized_geo/' + low_name + ".xyz")
    this_run.scrpath = (basic_path + 'scr/' + low_name +"/optim.xyz")
    this_run.outpath = (basic_path + 'outfiles/' + base_name + ".out")
    this_run.configure(metal,ox,spin,has_o2,lig1,lig2)
    this_run.ID = ID
    if os.path.exists(this_run.geopath):
        this_run.geo_exists = True
        this_run.obtain_mol3d()
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
                    print(lines)
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
                        this_run.obtain_mol3d()
                if not this_run.geo_exists:
                        try:
                                this_run.extract_geo()
                        except:
                                print("ERROR: scr not found for" + str(this_run.geopath))
        if not this_run.converged:
                print(' job  ' + str(this_run.outpath) + 'not converged')
                if not this_run.geo_exists:
                        try:
                                this_run.extract_geo()
                        except:
                                print("ERROR: scr not found for" + str(this_run.geopath))
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


def process_runs(bound_runs,unbound_runs):
    final_results=dict()
    matched = False
    number_of_matches  = 0
    print('processing all converged runs')
    for ID in unbound_runs.keys():
        matched = 0 
        unbound_run = unbound_runs[ID]
        this_name = unbound_run.name
        this_ID = ID
        if ID in bound_runs.keys():
            bound_run = bound_runs[ID]
            matched = True
            number_of_matches += 1
        if matched:
            print('matched ID: '+ str( ID) + ' files ' + str(bound_run.name) + ' and ' + str(unbound_run.name))
            final_results[this_ID] = Comp(this_ID)
            final_results[this_ID].ID = this_ID
            final_results[this_ID].boundenergy = str(float(bound_run.energy)*HF_to_Kcal_mol)
            final_results[this_ID].unboundenergy = str(float(unbound_run.energy)*HF_to_Kcal_mol)
            final_results[this_ID].lig1 = unbound_run.lig1
            final_results[this_ID].lig2 = unbound_run.lig2
            final_results[this_ID].boundss_act = bound_run.ss_act
            final_results[this_ID].unboundss_act = unbound_run.ss_act
            final_results[this_ID].boundss_target = bound_run.ss_target
            final_results[this_ID].unboundss_target = unbound_run.ss_target
        else:
            print('unmatched ID: '+ str( ID) + ' files ' + str(unbound_run.name)+ ' has no partner' )

    return final_results

full_list_of_props = ['name','converged','spin','lig1','lig2','bound','ss_act','ss_target','energy']
summary_list_of_props = ['ID','lig1','lig2','unboundenergy','boundenergy','unboundss_act','unboundss_target','boundss_act','boundss_target']
new_file = open('full_results.txt','w')
new_file.write(','.join(full_list_of_props) + '\n')
summary_file = open('summary_results.txt','w')
summary_file.write(','.join(summary_list_of_props) + '\n')
targetpaths=sorted(glob.glob("completejobs/*.done"))
unbound_runs = dict()
bound_runs = dict()
for files in targetpaths:
        this_run = test_terachem_go_convergence(files)
        if this_run.converged:
            extrct_props = scfextract(this_run,full_list_of_props)
            writeprops(extrct_props,new_file,1)
            if this_run.bound:
                bound_runs.update({this_run.ID:this_run})
            else:
                unbound_runs.update({this_run.ID:this_run})

targetpaths=sorted(glob.glob("jobs/*.in"))
for files in targetpaths:
        this_run = test_terachem_go_convergence(files)
        if this_run.converged:
            extrct_props = scfextract(this_run,full_list_of_props)
            writeprops(extrct_props,new_file,1)
            if this_run.bound:
                bound_runs.update({this_run.ID:this_run})
            else:
                unbound_runs.update({this_run.ID:this_run})

final_results = process_runs(bound_runs,unbound_runs)
for keys in sorted(final_results):
    this_result = final_results[keys]
    extrct_props = scfextract(this_result,summary_list_of_props)
    writeprops(extrct_props,summary_file,1)



new_file.close()
summary_file.close()
