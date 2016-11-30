import glob, string, os
from geometry import *
from atom3D import *
from globalvars import globalvars
from mol3D import*

HF_to_Kcal_mol = 627.503
EV_to_Kcal_mol = 23.06055

def name_breakdown(name):
    name_split = name.split('_')
    print(name_split)
    ret_dict = dict()
    ret_dict.update({'ID': int(name_split[0])})
    ret_dict.update({'metal':name_split[1]})
    ret_dict.update({'ox':int(name_split[2])})
    ret_dict.update({'eq_ind':int(name_split[4])})
    ret_dict.update({'ax1_ind':int(name_split[6])})
    ret_dict.update({'ax2_ind':int(name_split[8])})
    ret_dict.update({'spin':int(name_split[9])})
    return ret_dict
class mopac_run:
    """ This is a class for each run"""

    def __init__(self,name):

        self.name = name
        self.outpath  = 'undef'
        self.geopath = 'undef'
        self.geo_exists = False
        self.output_exists = False
        self.converged = 'N'
        self.time = 'undef'
        self.energy = 0
        self.spin = 'undef'
        self.eqlig_charge = 'undef'
        self.ssq = 0
        self.star = 0
        self.gene = 'undef'
    def obtain_mol3d(self,geopath,mol_geopath):
        this_mol = mol3D()
        this_mol.readfromxyz(geopath + '/' + self.name + '.xyz')
        this_mopac_mol = mol3D()
        this_mopac_mol.readfromxyz(mol_geopath + '/' + self.name + '.xyz')
        self.mol = this_mol
        self.momol = this_mopac_mol
        self.rmsd = this_mol.rmsd(this_mopac_mol)
#    def obtain_ML_dists(self):
    #    self.min_dist = minimum_ML_dist(self.mol)
    #    self.max_dist = maximum_ML_dist(self.mol)
    #def configure_ligands(self):
     #   this_ax_lig = liganddict[self.axlig]
     #   this_eq_lig = liganddict[self.eqlig]
     #   self.axlig_dent = this_ax_lig[0]
     #   self.eqlig_dent = this_eq_lig[0]
     #   self.axlig_charge = this_ax_lig[1]
      #  self.eqlig_charge = this_eq_lig[1]
      #  self.axlig_connect = this_ax_lig[2]
      #  self.eqlig_connect = this_eq_lig[2]
      #  self.axlig_natoms = this_ax_lig[3]
       # self.eqlig_natoms = this_eq_lig[3]
       # self.axlig_mdelen =  this_ax_lig[4]
       # self.eqlig_mdelen = this_eq_lig[4]
targetpaths=sorted(glob.glob("pm7test/outfiles/*.out"))
mol_geopath = "/home/jp/Dropbox/Main/optimal_mol_des/automated_search_one/keep/molpac"
geopath ="/home/jp/Dropbox/Main/optimal_mol_des/automated_search_one/keep/obs/pm7test/optimized_geo"
print("Extracing data from *.out files: " + str(len(targetpaths) )+  " file found\n"  )
print_list = dict()
spin_list = dict()
ID_list = dict()
rmsd_list = dict()
mp_energy = dict()
for targets in targetpaths:
    with open(targets) as f:
        conv_flag =  False
        this_name = (os.path.basename(targets)).strip('.out')
        this_ret_dict = name_breakdown(this_name)
        this_run = mopac_run(this_name)
        data=f.readlines()
        found_data = 0
        found_time = 0
        this_geo = list()
        for i,lines in enumerate(data):
            if str(lines).find('Optimization Converged.') != -1:
                conv_flag = True
            if str(lines).find('FINAL ENERGY') != -1:
                this_run.energy =str(float(lines.split()[2])*HF_to_Kcal_mol)
            if str(lines).find('SPIN S-SQUARED') != -1:
                this_str=(lines.split())
                this_run.ssq =float( this_str[2])
                this_run.star = float(this_str[4].strip('()'))
        if conv_flag: 
            print(this_name + ' converged, checking geo')
            this_run.obtain_mol3d(geopath,mol_geopath)
            print_list.update({this_name:this_run.energy})
            spin_list.update({this_name:this_ret_dict['spin']})
            ID_list.update({this_name:this_ret_dict['ID']})
            rmsd_list.update({this_name:this_run.rmsd})
            with open("/home/jp/Dropbox/Main/optimal_mol_des/automated_search_one/keep/molpac/" + this_name + '.out') as mf:
                mdata=mf.readlines()
                for j,mlines in enumerate(mdata):
                    if str(mlines).find('TOTAL ENERGY') != -1:
                        #:wsprint(lines)
                        this_run.menergy =str(float(mlines.split()[3])*EV_to_Kcal_mol)
                        print('found mopac for ' + this_name)
                        mp_energy.update({this_name:this_run.menergy})
with open('dft_results.csv','w') as f:
    for keys in print_list.keys():
        stw =str(ID_list[keys]) + ','+ str(spin_list[keys]) + ','+ str(keys) + ','+ str(print_list[keys])+','+ str(rmsd_list[keys]) + ', ' + str(mp_energy[keys])+'\n'
        f.write(stw)
# at this stage, all of the real data has b



