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
HF_to_Kcal_mol = 627.509###
###########################

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



class DFTRun:
    """ This is a class for each run"""
    numRuns = 0
    def __init__(self,name):
        self.numRuns += 1
        self.name = name
        self.outpath  = 'undef'
        self.geopath = 'undef'
        self.init_geopath = 'undef'
        self.progpath  = 'undef'
        self.geo_exists = False
        self.output_exists = False
        self.converged = 'N'
        self.time = 'undef'
        self.energy = 0
        self.spin = 'undef'
        self.eqlig_ind = 'undef'
        self.axlig1_ind = 'undef'
        self.axlig2_ind = 'undef'
        self.eqlig = 'undef'
        self.axlig = 'undef'
        self.axlig = 'undef'

        self.ss_target = 0
        self.ss_act = 0
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
    def extract_prog(self):
            cmd_str = 'python optgeo_extract.py '+ self.scrpath + ' ' + self.progpath
 #           print(cmd_str)
            subprocess.call(cmd_str,shell=True)
    def extract_geo(self):
            cmd_str = 'python optgeo_extract.py '+ self.scrpath + ' ' + self.geopath
  #          print(cmd_str)
            subprocess.call(cmd_str,shell=True)

    def obtain_ML_dists(self):
        self.min_dist = minimum_ML_dist(self.mol)
        self.max_dist = maximum_ML_dist(self.mol)
    def configure(self,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat):
        self.metal = metal
        self.ox = ox
        self.spin = spin
        self.eqlig_ind = eqlig
        self.axlig1_ind = axlig2
        self.axlig2_ind = axlig1
        ligands_dict = get_ligands()
        self.eqlig = ligands_dict[int(self.eqlig_ind)][0]
        self.axlig1 = ligands_dict[int(self.axlig1_ind)][0]
        self.axlig2 = ligands_dict[int(self.axlig2_ind)][0]

class Comp:
    """ This is a class for each unique composition and configuration"""
    def __init__(self,name):
        self.name = name
        self.ox =' undef'
        self.metal= 'undef'
        self.axlig1 = 'undef'
        self.axlig2 = 'undef'
        self.eqlig = 'undef'
        self.axlig1_ind = 'undef'
        self.axlig2_ind = 'undef'
        self.eqlig_ind = 'undef'
        self.alpha = 'undef'
        self.HSenergy = 'undef'
        self.LSenergy= 'undef'
        self.times = list()
        self.splitenergy = 0

    def process(self):
        self.splitenergy = str((float(self.HSenergy) - float(self.LSenergy))*HF_to_Kcal_mol)
    def find_fitness(self):
        ref_value = 15.0
#        print(self.splitenergy)
        en =-1*numpy.power((float(self.splitenergy)/ref_value),2.0)
#        print(en)
        self.fitness = numpy.exp(en)

def writeprops(extrct_props,startpoints,newfile,do_strip):
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
    path_dictionary = setup_paths()
    basic_path = get_run_dir()
    ID,low_name,base_name,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat,gene = translate_job_name(job)
    ### flag
    converged =  False
    geo_exists = False
    ### test if geo exits
    this_run=DFTRun(base_name)
    this_run.converged = False
    this_run.geo_exists = False
    this_run.geopath = (path_dictionary["optimial_geo_path" ] + base_name + ".xyz")
    this_run.progpath = (path_dictionary["progress_geo_path" ] + base_name + ".xyz")

    this_run.outpath = (path_dictionary["geo_out_path" ] + base_name + ".out")
    this_run.scrpath = (basic_path + 'scr/geo/' + base_name +"/optim.xyz")
    this_run.inpath = (basic_path + 'jobs/' + base_name +".in")
    this_run.comppath = (basic_path + 'completejobs/' + base_name +".in")

    this_run.configure(metal,ox,eqlig,axlig1,axlig2,spin,spin_cat)
    this_run.gene = gene
    this_run.ID = ID
    this_run.spin_cat = spin_cat
    print('run is set up')
    if os.path.exists(this_run.geopath):
        this_run.geo_exists = True
    if os.path.exists(this_run.outpath):
        ### file is found, check if converged
        with open(this_run.outpath) as f:
            data=f.readlines()
            found_conv =False 
            found_data =False
            found_time = False 
            for i,lines in enumerate(data):
                if str(lines).find('Optimization Converged.') != -1:
                   found_conv = True
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
                print('run converged ' + str(this_run.name) + ' and now testing geoex ' )
                if this_run.geo_exists:
                        print('geo exists for ' + this_run.name)
                if not this_run.geo_exists:
                        print('geoex not found at ' +str(this_run.geopath) +  ' for ' + this_run.name)
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
                if not os.path.exists(this_run.comppath):
                        print('this run does not have finished filese')
                        shutil.copy(this_run.inpath,this_run.comppath)
                        logger(path_dictionary['state_path'],str(datetime.datetime.now()) + " moving  " + str(this_run.name) + " to " + str(this_run.comppath))
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
    return this_run
def test_terachem_sp_convergence(job):
    ### get paths
    path_dictionary = setup_paths()
    ## get job properties
    base_name = os.path.basename(job)
    base_name = base_name.strip('.in')
    ID,low_name,base_name,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat,gene = translate_job_name(job)

    ### flag
    converged =  False
    ### test if geo exits
    this_run=DFTRun(base_name)
    this_run.outpath = (path_dictionary["vertical_out_path" ]  + base_name + ".out")
    print("checking ",this_run.outpath)
    this_run.configure(metal,ox,eqlig,axlig1,axlig2,spin,spin_cat)
    this_run.gene = gene
    this_run.ID = ID
    this_run.spin_cat = spin_cat

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
                    this_run.ss_act =float( this_str[2])
                    this_run.ss_target = float(this_str[4].strip('()'))
        if (found_data == True) and (found_time == True) and (found_conv == True):
            this_run.converged = True
    return this_run

def process_runs_sp(LS_runs,HS_runs):
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
            final_results[this_ID].LSenergy = str(float(LS_run.energy))
            final_results[this_ID].HSenergy = str(float(HS_run.energy))
            final_results[this_ID].process()
            final_results[this_ID].eqlig_ind = LS_run.eqlig_ind
            final_results[this_ID].axlig1_ind = LS_run.axlig1_ind
            final_results[this_ID].axlig2_ind = LS_run.axlig2_ind
            final_results[this_ID].eqlig_ind = LS_run.eqlig
            final_results[this_ID].axlig1_ind = LS_run.axlig1
            final_results[this_ID].axlig2_ind = LS_run.axlig2
            final_results[this_ID].HSss_act = HS_run.ss_act
            final_results[this_ID].LSss_act = LS_run.ss_act
            final_results[this_ID].LSss_target = LS_run.ss_target
            final_results[this_ID].HSss_target = HS_run.ss_target
        else:
            print('unmatched ID: '+ str( ID) + ' files ' + str(LS_run.name)+ ' has no partner' )
    return final_results

def process_runs_geo(LS_runs,HS_runs):
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
            final_results[this_ID].LSenergy = str(float(LS_run.energy))
            final_results[this_ID].HSenergy = str(float(HS_run.energy))
            final_results[this_ID].process()
            final_results[this_ID].eqlig_ind = LS_run.eqlig_ind
            final_results[this_ID].axlig1_ind = LS_run.axlig1_ind
            final_results[this_ID].axlig2_ind = LS_run.axlig2_ind
            final_results[this_ID].eqlig_ind = LS_run.eqlig
            final_results[this_ID].axlig1_ind = LS_run.axlig1
            final_results[this_ID].axlig2_ind = LS_run.axlig2
            final_results[this_ID].HSss_act = HS_run.ss_act
            final_results[this_ID].LSss_act = LS_run.ss_act
            final_results[this_ID].LSss_target = LS_run.ss_target
            final_results[this_ID].HSss_target = HS_run.ss_target
            final_results[this_ID].LS_rmsd = LS_run.rmsd
            final_results[this_ID].HS_rmsd = HS_run.rmsd
        else:
            print('unmatched ID: '+ str( ID) + ' files ' + str(LS_run.name)+ ' has no partner' )
    return final_results






