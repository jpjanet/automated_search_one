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
import shutil
from geometry import *
from atom3D import *
from globalvars import globalvars
from mol3D import*



########### UNIT CONVERSION
HF_to_Kcal_mol = 627.503

ligandsforexp = [['doxy1',['lib','O2','O',-1]],
                 ['thiocyanate',['lib','SCN','S',-1]],
                 ['Cl',['halo','Cl','Cl',-1]],
                 ['hydroxyl',['lib','OH','O',-1]],
                 ['water',['lib','H2O','O',0]],
                 ['isothiocyanate',['lib','NCS','N',-1]],
                 ['ammonia',['lib','NH3','N',0]],
                 ['cyanide',['lib','CN','C',-1]],
                 ['carbonyl',['lib','CO','C',-1]]]
liganddict = {'scn':[1,-1,'S',3,0.03],'cl':[1,-1,'Cl',1,0],'c2h3ns':[1,0,'C',6,-0.49],'h2o':[1,0,'O',3,1.24],
              'ncs':[1,-1,'N',3,0.49],'nh3':[1,0,'N',4,0.84],'cn':[1,-1,'C',2,-0.49],'co':[1,0,'C',2,-0.89],'bipy':[2,0,'N',20,0.49],
              'phen':[2,0,'N',22,0.49],'ox':[2,-2,'O',6, 0.89],'acac':[2,0,'O',15,0.89],'en':[2,0,'N',12,0.84],
              'porphyrin':[4,-2,'N',36,0.49],'pisc':[1,0,'C',25,-0.49],'tbuc':[2,-2,'O',24,0.89]}

metal_ox_states = [['Cr',[2,3]],['Mn',[2,3]],['Fe',[2,3]],['Co',[2,3]],['Ni',[2]]]

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



class Run:
    """ This is a class for each run"""
    numRuns = 0
    def __init__(self,name):
        Run.numRuns += 1
        self.name = name
        self.path  = 'undef'
        self.converged = 'N'
        self.time = 'undef'
        self.energy = 0
        self.spin = 'undef'
        self.ox =' undef'
        self.metal= 'undef'
        self.axlig = 'undef'
        self.eqlig = 'undef'
        self.alpha = 'undef'
        self.runnum = 'undef'
        self.time = 'undef'
        self.spin_cat = 'undef'
        self.axlig_connect = 'undef'
        self.eqlig_connect = 'undef'
        self.axlig_natoms = 'undef'
        self.eqlig_natoms = 'undef'
        self.axlig_dent = 'undef'
        self.eqlig_dent = 'undef'
        self.axlig_mdelen = 'undef'
        self.eqlig_mdelen = 'undef'
        self.axlig_charge = 'undef'
        self.eqlig_charge = 'undef'
        self.ssq = 0
        self.star = 0
        self.gtx970 = 0

    def obtain_mol3d(self,geopath):
        this_mol = mol3D()
        this_mol.readfromxyz(geopath + self.name + '.xyz')
        self.mol = this_mol
    def obtain_ML_dists(self):
        self.min_dist = minimum_ML_dist(self.mol)
        self.max_dist = maximum_ML_dist(self.mol)
    def configure_ligands(self):
        this_ax_lig = liganddict[self.axlig]
        this_eq_lig = liganddict[self.eqlig]
        self.axlig_dent = this_ax_lig[0]
        self.eqlig_dent = this_eq_lig[0]
        self.axlig_charge = this_ax_lig[1]
        self.eqlig_charge = this_eq_lig[1]
        self.axlig_connect = this_ax_lig[2]
        self.eqlig_connect = this_eq_lig[2]
        self.axlig_natoms = this_ax_lig[3]
        self.eqlig_natoms = this_eq_lig[3]
        self.axlig_mdelen =  this_ax_lig[4]
        self.eqlig_mdelen = this_eq_lig[4]
#        print(this_ax_lig)
class Comp:
    """ This is a class for each unique composition and configuration"""
    def __init__(self,name):
        self.name = name
        self.ox =' undef'
        self.metal= 'undef'
        self.axlig = 'undef'
        self.eqlig = 'undef'
        self.alpha = 'undef'
        self.HSenergy = 'undef'
        self.LSenergy= 'undef'
        self.runnumbers  = list()
        self.times = list()
        self.axlig_connect = 'undef'
        self.eqlig_connect = 'undef'
        self.axlig_natoms = 'undef'
        self.eqlig_natoms = 'undef'
        self.axlig_dent = 'undef'
        self.eqlig_dent = 'undef'
        self.axlig_mdelen = 'undef'
        self.eqlig_mdelen = 'undef'
        self.axlig_charge = 'undef'
        self.eqlig_charge = 'undef'
        self.hs_max_dist = 'undef'
        self.hs_min_dist = 'undef'
        self.ls_max_dist = 'undef'
        self.ls_min_dist = 'undef'
        self.splitenergy = 0

    def process(self):
        self.splitenergy = str((float(self.HSenergy) - float(self.LSenergy))*HF_to_Kcal_mol)

def dash_advance(string,break_ind):
    break_ind =string.find('_',break_ind+1)
    if break_ind == -1:
##        print("file string not formatted correctly")
##        print(string)
        return 1
    else:
        return break_ind



def stringteardown(string):
    """   :rtype: Run   """
    shortname = os.path.splitext(os.path.basename(string))[0]
    thisrun = Run(shortname)
    thisrun.path = string
    break_ind =0
##    print('extracting ' + shortname)
    break_ind = dash_advance(shortname,break_ind)
    old_break = break_ind
    thisrun.runnum =int(shortname[:break_ind])
    break_ind = dash_advance(shortname,old_break)
    thisrun.metal = shortname[(old_break + 1):break_ind]
    old_break = break_ind

    break_ind = dash_advance(shortname,old_break)
    thisrun.ox = shortname[(old_break + 1):break_ind]
    old_break = break_ind

    break_ind = dash_advance(shortname,old_break)
    thisrun.spin = shortname[(old_break + 1):break_ind]
    old_break = break_ind

    break_ind = dash_advance(shortname,old_break)
    old_break = break_ind
    break_ind = dash_advance(shortname,old_break)
    next_ind = dash_advance(shortname,break_ind)
    thisrun.axlig = shortname[(old_break + 1):(break_ind)]
    old_break = break_ind


    break_ind = dash_advance(shortname,old_break)
    old_break = break_ind
    break_ind = dash_advance(shortname,old_break)
    next_ind = dash_advance(shortname,break_ind)
    thisrun.eqlig = shortname[(old_break + 1):(break_ind)]
    old_break = break_ind


    thisrun.alpha = shortname[(break_ind+1):(break_ind+4)]
    return thisrun

def writeprops(extrct_props,startpoints,newfile,do_strip):
#    for word in xrange(0,len(extrct_props)):
 #       if do_strip == 0:
 #           wbuffer=str(extrct_props[word]).strip()
 #       else:
 #           wbuffer=str(extrct_props[word])
 #       wbuffer=wbuffer.ljust(startpoints[word + 1])
#        print(wbuffer)
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
fi =dict()
pen = 0
list_of_props = ['name','runnum','converged','spin','ox','axlig','eqlig','axlig_charge','eqlig_charge',
                 'axlig_dent','eqlig_dent','axlig_connect','eqlig_connect','axlig_natoms','eqlig_natoms','axlig_mdelen','eqlig_mdelen'
                 ,'alpha','min_dist','max_dist','energy','gtx970','time','ssq','star']
summary_list_of_props = ['name','metal','ox','alpha','rmsd','axlig','eqlig','axlig_charge','eqlig_charge',
                 'axlig_dent','eqlig_dent','axlig_connect','eqlig_connect','axlig_natoms','eqlig_natoms','axlig_mdelen','eqlig_mdelen'
                 ,'hs_min_dist','hs_max_dist','ls_min_dist','ls_max_dist','splitenergy','ssq','star']

newfile = open('cf_full.txt','w')
#newfile.write('Name        NUM            CONV  SPIN  OX     TIME ALPHA     MINDIST        MAXDIST    ENERGY\n')
newfile.write(','.join(list_of_props) + '\n')
summary_file = open('data_summary.txt','w')
#summary_file.write("NAME         METAL             AXLIG         EQLIG               OX                 ALPHA              RMSD              DELTA(HS,LS)\n")
summary_file.write(','.join(summary_list_of_props) + '\n')

#                  123456789012345
#                                   1234567890
#                                             12345
#                                                  12345678901
#                                                             1234567890
#                                                                       1234567890123456
#                                                                                       1234567890
# 6789012345678901234567890123456789012345678901234567890132456798012354678901235467890
summary_file_MR = open('data_summary_MR.dat','w')
column_labels = ['NAME','METAL','LIG','OX','ALPHA','DELTA(HS,LS)']
row_label_file = open('row_labels.txt','w')
column_label_file = open('column_labels.txt','w')
for labels in column_labels:
    column_label_file.write(labels + "\n")
column_label_file.close()
startpoints = [0,40,6,6,6,6,6,6,6,6,6,6,15,20,20,20,20,20,20,10] # These are actually lenghts. These should be made dynamic based on the properties requested
extrct_props = [1,2,3,4,5]
summary_startpoints1 = [0,30,30,30,30,30,306,6,6,20,20,20,20,20,20,20,20,20,20,20,20,20]
summary_startpoints_MR = [10]*13
newnames=dict()
shortname=[]
all_runs=dict();
targetpaths=sorted(glob.glob("outfiles/*.out"))
geopath = "optimized_geo/"
print("Extracing data from *.out files: " + str(len(targetpaths) )+  " file found\n"  )

for j,paths in enumerate(targetpaths):
    this_run=stringteardown(paths)
    with open(this_run.path) as f:
        data=f.readlines()
        found_data = 0
        found_time = 0
        unconv = 1 #false
        for i,lines in enumerate(data):
            if str(lines).find('  Device 0:      GeForce GTX 970') != -1:
#                print('GTX')
#                print(lines)
                this_run.gtx970 = 1
            if str(lines).find('Converged!') != -1:
                unconv = 0
#                print('conv!')
            if str(lines).find('Optimization Converged.') != -1:
                unconv = 0
#                print('conv!')

            if str(lines).find('FINAL ENERGY') != -1:
#                print(lines)
                this_run.energy =str(lines.split()[2])
#                print(this_run.energy)
                found_data = 1
            if str(lines).find('Total processing time') != -1:
                this_run.time=str(lines.split()[3])
#                print(lines)
                found_time =1
            if str(lines).find('SPIN S-SQUARED') != -1:
                this_str=(lines.split())
                this_run.ssq =float( this_str[2])
                this_run.star = float(this_str[4].strip('()'))

    if (found_data == 1) and (found_time == 1) and (unconv == 0):
        this_run.converged = 'Y'
#        print('run converged')
    try:
        this_run.obtain_mol3d(geopath)
        this_run.obtain_ML_dists()
    except:
        this_run.min_dist = 0
        this_run.max_dist = 0

    this_run.configure_ligands()
    extrct_props = scfextract(this_run,list_of_props)
 #   print(extrct_props)
    writeprops(extrct_props,startpoints,newfile,1)
    all_runs[str(this_run.runnum)+'_'+str(this_run.alpha)]  = this_run
# at this stage, all of the real data has been imported. Now, runs must be totaled, averaged and written

print("All files read! # of files = " + str(len(all_runs)) +" A full set of results is avaiable in SCF_full.txt")
print("Results will now be post-processed, added and compared. Stand by...")
final_results=dict()
matched = 0
old_metal = " "
old_lig = " "
old_metal = " "
old_ox = " "
unproc_runs = all_runs
number_of_matches  = 0
for run_names in unproc_runs.keys():
    matched = 0 
#    print('\n')
 #   print(unproc_runs[run_names].name)
    if (unproc_runs[run_names].runnum % 2):
       # print('Odd numbered run, one of the pair')
           # find matching partner
        if (unproc_runs[run_names].runnum < 900):
            try:
                partner_name = str(unproc_runs[run_names].runnum+1) + '_' + unproc_runs[run_names].alpha
     #          print(run_names + ' and ' + partner_name)
                HSpartner = all_runs[partner_name]
                LSpartner = all_runs[run_names]
                matched = 1
            except:
                print('run ' + str(run_names) +  ' appears unmatched')
        if (unproc_runs[run_names].runnum >= 900):
            try:
                partner_name = str(unproc_runs[run_names].runnum-1) + '_' + unproc_runs[run_names].alpha
     #          print(run_names + ' and ' + partner_name)
                LSpartner = all_runs[partner_name]
                HSpartner = all_runs[run_names]
                matched = 1
            except:
                print('run ' + str(run_names) +  ' appears unmatched')

        if ((matched) and (LSpartner.converged == "Y") and (HSpartner.converged == "Y")  and (LSpartner.energy != 0) and (HSpartner.energy != 0)):
                name = (unproc_runs[run_names].metal + '('  +  str(unproc_runs[run_names].ox) + ')'  + '_ax_' + unproc_runs[run_names].axlig
                                                    + '_eq_'  + unproc_runs[run_names].eqlig + str(unproc_runs[run_names].alpha))

                final_results[name] = Comp(name)
                numlist = [unproc_runs[run_names].runnum,unproc_runs[run_names].runnum +1 ],
                number_of_matches += 1
                if not (LSpartner.metal ==  HSpartner.metal) or  not (LSpartner.axlig == HSpartner.axlig) or not (LSpartner.alpha == HSpartner.alpha):
                    print('error, partner mismatch between ' + str(numlist))
                if (HSpartner.spin < LSpartner.spin):
                    print(' possible spin mismatch at  ')
                    print(str(run_names))
                final_results[name].metal = unproc_runs[run_names].metal
                final_results[name].ox = unproc_runs[run_names].ox
                final_results[name].axlig = unproc_runs[run_names].axlig
                final_results[name].eqlig = unproc_runs[run_names].eqlig
                final_results[name].axlig_dent = unproc_runs[run_names].axlig_dent
                final_results[name].eqlig_dent = unproc_runs[run_names].eqlig_dent
                final_results[name].axlig_charge = unproc_runs[run_names].axlig_charge
                final_results[name].eqlig_charge = unproc_runs[run_names].eqlig_charge
                final_results[name].axlig_natoms = unproc_runs[run_names].axlig_natoms
                final_results[name].eqlig_natoms = unproc_runs[run_names].eqlig_natoms
                final_results[name].axlig_connect = unproc_runs[run_names].axlig_connect
                final_results[name].eqlig_connect = unproc_runs[run_names].eqlig_connect
                final_results[name].axlig_mdelen = unproc_runs[run_names].axlig_mdelen
                final_results[name].eqlig_mdelen = unproc_runs[run_names].eqlig_mdelen
                final_results[name].runnumbers = numlist
                final_results[name].alpha = unproc_runs[run_names].alpha
                final_results[name].HSenergy = HSpartner.energy
                final_results[name].hs_max_dist = HSpartner.max_dist
                final_results[name].hs_min_dist = HSpartner.min_dist
                final_results[name].ls_max_dist = LSpartner.max_dist
                final_results[name].ls_min_dist = LSpartner.min_dist
                final_results[name].LSenergy = LSpartner.energy
                try:
                    shutil.copy(geopath + LSpartner.name +'.xyz','finalgeos/')
                    shutil.copy(geopath + HSpartner.name +'.xyz','finalgeos/')
                except:
                    print(LSpartner.name,'not found geo')

                LS_ss_error = abs(LSpartner.star - LSpartner.ssq)
                HS_ss_error = abs(HSpartner.star - HSpartner.ssq)

                try:
                    final_results[name].rmsd = HSpartner.mol.rmsd(LSpartner.mol)
                    final_results[name].process()
                except:
                    final_results[name].rmsd = 0
                    final_results[name].splitenergy = 0

 #           print(unproc_runs[run_names].name)
                if (LS_ss_error >= HS_ss_error):
                    final_results[name].star = LSpartner.star
                    final_results[name].ssq = LSpartner.ssq
                else:
                    final_results[name].star = HSpartner.star
                    final_results[name].ssq = HSpartner.ssq

        #    if (unproc_runs[run_names].runnum == 61):
        #        print(name)
        #        print('hse is ' + str(final_results[name].HSenergy))
        #        print('lse is ' + str(final_results[name].LSenergy))
        #        print('delta is ' + str(final_results[name].splitenergy))


print("\n\n")
indiv_files = dict()
for m_ox in metal_ox_states:
    metal =  m_ox[0]
    for ox in m_ox[1]:
        name = metal[0].lower() + str(ox) + 'f'
row_labels = []
for keys in sorted(final_results):
    this_result = final_results[keys]
    extrct_props = scfextract(this_result,summary_list_of_props)
    writeprops(extrct_props,summary_startpoints1,summary_file,1)
  # Now construct machine-readable file and labels
  #  writeprops(extrct_props,summary_startpoints_MR,summary_file_MR,0)
#for labels in row_labels:
 #   row_label_file.write(labels + "\n")
#row_label_file.close()

newfile.close()

summary_file.close()
#summary_file_MR.close()


    
    
print("All work complete. The summary results can be found in SCF_summary.txt \n")
print("A machine-readable file is provided as SCF_summary_MR.dat, with labels in row_labels.txt and column_labels.txt")


