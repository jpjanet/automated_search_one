import glob, string, os
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
    #def obtain_mol3d(self,geopath):
    #    this_mol = mol3D()
    #    this_mol.readfromxyz(geopath + self.name + '.xyz')
    #    self.mol = this_mol
    #def obtain_ML_dists(self):
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
targetpaths=sorted(glob.glob("keep/molpac/*.out"))
geopath = "mopac_geooptimized_geo/"
print("Extracing data from *.out files: " + str(len(targetpaths) )+  " file found\n"  )
print_list = dict()
spin_list = dict()
ID_list = dict()
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
        in_cord = False
        for i,lines in enumerate(data):
            if str(lines).find('TOTAL ENERGY') != -1:
                #:wsprint(lines)
                this_run.energy =str(float(lines.split()[3])*EV_to_Kcal_mol)
            if str(lines).find('Converged!') != -1:
                unconv = 0
#                print('conv!')
            if str(lines).find('SCF FIELD WAS ACHIEVED ') != -1:
                conv_flag =  True
                this_run.converged = True
#                print('conv!')

            if (str(lines).find('CARTESIAN COORDINATES') != -1) and  (conv_flag):
                in_cord = True
                print('found final geo')

            if in_cord:
                if (str(lines).find('Empirical Formula') != -1):
                    in_cord =  False
                    print('end of geo')
                else:
                    if lines.strip():
                        this_geo.append(list(lines[1:]))
#        print(this_geo)
        with open('keep/molpac/'+this_name+'.xyz','w') as f:
            f.write(str(int(len(this_geo))-1)+'\n')
            f.write('#'+this_name+'\n')
            for i,elements in enumerate(this_geo):
                if not i==0:
                    line_tw = ''.join(elements[5:])
                    line_tw= line_tw.lstrip()
                    f.write(line_tw)
        print_list.update({this_name:this_run.energy})
        spin_list.update({this_name:this_ret_dict['spin']})
        ID_list.update({this_name:this_ret_dict['ID']})
with open('mopac_results.csv','w') as f:
    for keys in print_list.keys():
        stw =str(ID_list[keys]) + ','+ str(spin_list[keys]) + ','+ str(keys) + ','+ str(print_list[keys])+'\n'
        f.write(stw)
# at this stage, all of the real data has b



