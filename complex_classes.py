import glob
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil


from prep_calc import * 

class TM_complex:
    def __init__(self,ID,core,oxidation_state):
        self.free_sites = [1,2,3,4,5,6]
        ### mark 3x bidentate
        self.three_bidentate = False
        self.ID = ID
        self.geo = 'oct'
        self.ready_for_assembly = False
        self.ligands_dict=get_ligands()
        self.core  = core
        self.oxidation_state = int(oxidation_state)
        self.ax_dent =  False
        self.eq_dent = False
        self.eq_ligands = list()
        self.ax_ligands= list()
        self.ax_inds = list()
        self.eq_inds = list()
    def random_gen(self):
        ## set equitorial first
        self._get_random_equitorial()
        self._get_random_axial()
        self._name_self()
    def check_self(self):
        if 0 in self.ax_inds:
            self.geo ='spy'
        else:
            self.geo = 'oct'
    def check_self_loud(self):
#        print('checking in ' + str(self.ax_inds)  )
        if 0 in self.ax_inds:
            self.geo ='spy'
            print('spy')
        else:
            self.geo = 'oct'
  #          print('oct')


    def _name_self(self):
        self.check_self()
        self.name =(str(self.core) + '_' + str(self.oxidation_state)+'_eq_' + str(self.eq_inds[0]) + '_ax1_' + str(self.ax_inds[0]) + '_ax2_' + str(self.ax_inds[1])) 

    def copy(self,partner):
         self.ax_dent = partner.ax_dent
         self.ax_inds = partner.ax_inds
         self.ax_ligands = partner.ax_ligands
         self.eq_dent = partner.eq_dent
         self.eq_inds = partner.eq_inds
         self.eq_oc = partner.eq_oc
         self.eq_ligands = partner.eq_ligands
         self.three_bidentate = partner.three_bidentate
         self._name_self()
    def _get_random_equitorial(self):
        ### choose the equitorial ligand
        n = len(self.ligands_dict)
        eq_ind = 0
        while (eq_ind == 0):
            # eq site cannot be empty
            eq_ind = numpy.random.randint(low = 0,high = n)
            if (self.ligands_dict[eq_ind][0] == 'dppe') or (self.ligands_dict[eq_ind][0] == 'pph3'):
                eq_ind = 0
                print('preventing dppe/pph3 in equitorial')

        eq_ligand_properties  = self.ligands_dict[eq_ind][1]
        self.eq_dent = eq_ligand_properties[0]
        self.eq_oc  = int(4/self.eq_dent)
        self.eq_ligands = [self.ligands_dict[eq_ind][0] for i in range(0,self.eq_oc)]
        self.eq_inds = [eq_ind]
    def  _get_random_axial(self):
        ### choose axial ligands:
        n = len(self.ligands_dict)
        self.ready_for_assembly =  False
        self.ax_ligands = list()
        self.ax_inds = list()
        has_zero = False
        while not self.ready_for_assembly:
            ax_ind = numpy.random.randint(low = 0,high = n)
            #print(ax_ind)
            ax_ligand_properties  = self.ligands_dict[ax_ind][1]
            ax_dent = ax_ligand_properties[0]
            if ax_dent > 1 and not has_zero:
                if ((self.eq_dent == 2) and (ax_dent == 2) and (len(self.ax_ligands) == 0)):
                    three_bidentate  = True
                    self.ax_ligands = [self.ligands_dict[ax_ind][0],self.ligands_dict[ax_ind][0]]
                    self.ax_inds = [ax_ind, ax_ind]
                    self.ready_for_assembly = True
                    self.ax_dent = 2
                    self.ax_oc = [2]
                else:
                    self.ready_for_assembly = False
            elif ax_dent ==1:
                if ((ax_ind == 0) and (not has_zero)):
                    has_zero = True
                    self.ax_ligands.append(self.ligands_dict[ax_ind][0])
                    self.ax_inds.append(ax_ind)
                    self.geo = "spy"
                elif ax_ind !=  0:
                    self.ax_ligands.append(self.ligands_dict[ax_ind][0])
                    self.ax_inds.append(ax_ind)
                if (len(self.ax_ligands) ==2):
                    self.ax_dent = 1
                    self.ax_oc = [1,1]
                    self.ready_for_assembly = True
                else:
                     self.ready_for_assembly = False
        if self.geo == "spy":
            non_zero_lig = [i for i in self.ax_inds if (i != 0)]
            non_zero_lig.append(0)
            #print('nonzero ligs  = ',non_zero_lig)
            self.ax_inds  = non_zero_lig
            #print(self.ax_inds)
            ligs =[self.ligands_dict[j][0] for j in self.ax_inds]
            #print(ligs)
            self.ax_ligands = ligs
        self.ax_inds = sorted(self.ax_inds)
        self.ax_ligands = list()
        for ind in self.ax_inds:
            self.ax_ligands.append(self.ligands_dict[ind][0])
        self._name_self()

    def examine(self):
        print("name is " + self.name)
        print("eq", self.eq_ligands, self.eq_inds)
        print("axial", self.ax_ligands, self.ax_inds)

    def encode(self,ll):
        self.random_gen()
        ll = [int(item) for item in ll]
        self.replace_equitorial([ll[0]])
        self.replace_axial(ll[1:])
        self._name_self()

    def replace_equitorial(self,new_eq_ind):
        #prin?t('in repcoding, setting eq to ' + str(new_eq_ind))
        eq_ligand_properties  = self.ligands_dict[new_eq_ind[0]][1]
        self.eq_dent = eq_ligand_properties[0]
        self.eq_oc  = int(4/self.eq_dent)
        self.eq_ligands = [self.ligands_dict[new_eq_ind[0]][0] for i in range(0,self.eq_oc)]
        self.eq_inds = new_eq_ind
        if (self.ax_dent == 1) or ((self.ax_dent == 2) and (self.eq_dent ==2)):
                ## everything is ok!
                if (self.ax_dent == 2):
                    self.three_bidentate = True
        else: ## this complex cannot exist. keeping  equitorial,
              ## regenerating axial ligands
           print("complex with" + str(self.eq_ligands) + " and " + str(self.ax_ligands) + " cannot exist, randomizing axial")
           self._get_random_axial()
        self._name_self()

    def replace_axial(self,new_ax_ind):
        self.ax_ligands = list()
        self.ax_inds = list()
        n = len(self.ligands_dict)
        if 0 in new_ax_ind:
            self.geo ="spy"
        #print('new_ax_ind:'+ str(new_ax_ind))
        self.three_bidentate =  False
        for i,indices in enumerate(new_ax_ind):
            #print('trying to add ',indices )
            ax_ligand_properties = self.ligands_dict[indices][1]
            ax_dent = ax_ligand_properties[0]
            #print('ax dent ' + str(ax_dent))
            if (ax_dent > 1):
                if (ax_dent == 2) and (i == 0):
                    three_bidentate  = True
                    self.ax_ligands = [self.ligands_dict[indices][0],self.ligands_dict[indices][0]]
                    self.ax_inds = [indices, indices]
                    self.ax_dent = 2
                    self.three_bidentate = True
                    break
                else:
                    print('impossible, ax_dent  = ' + str(ax_dent))
                    break
            else:
                self.ax_ligands.append(self.ligands_dict[indices][0])
                self.ax_inds.append(indices)
                self.ax_dent = 1

        if (self.three_bidentate) and not (self.eq_dent == 2):
            ## this complex cannot exist
            ## regenerating equitorial ligands
            print("complex with" + str(self.eq_ligands) + " and " + str(self.ax_ligands) +
                  " cannot exist, regenerating equitorial ligands")
            self.ready_for_assembly  = False
            while not self.ready_for_assembly:
                eq_ind = numpy.random.randint(low = 0,high = n)
                eq_ligand_properties  = self.ligands_dict[eq_ind][1]
                eq_dent = eq_ligand_properties[0]
                if (eq_dent == 2):
                    self.eq_dent = 2
                    self.eq_inds = [eq_ind]
                    self.eq_oc = int(4/eq_dent)
                    self.eq_ligands = [self.ligands_dict[eq_ind][0] for i in range(0,self.eq_oc)]
                    self.ready_for_assembly =  True
        self.check_self_loud()
        if self.geo == "spy":
            #print(self.ax_inds)
            non_zero_lig = [i for i in self.ax_inds if (i != 0)]
            non_zero_lig.append(0)
            #print(non_zero_lig)
            self.ax_inds  = non_zero_lig
            #print(self.ax_inds)
            ligs =[self.ligands_dict[j][0] for j in self.ax_inds]
            #print(ligs)
        self.ax_inds = sorted(self.ax_inds)
        self.ax_ligands = list()
        for ind in self.ax_inds:
            self.ax_ligands.append(self.ligands_dict[ind][0])

        self._name_self()
        #print('final_ax_ind:'+ str(self.ax_inds))

    def exchange_ligands(self,partner,eq_swap):
        child = octahedral_complex(self.ligands_dict)
        child.copy(self) # copies this parent
        print("swapping from",partner.name," to ",self.name)
        self.examine()
        if eq_swap:
            print("swapping equitorial " + str(child.eq_inds) + ' -> ' + str(partner.eq_inds))
            child.replace_equitorial(partner.eq_inds)
        else:
            print("swapping axial"+ str(child.ax_inds) + ' -> ' + str(partner.ax_inds))
            child.replace_axial(partner.ax_inds)
        child.examine()
        child._name_self()
        return child


    def mutate(self):
        ## mutates either the axial
        ## or equitorial ligand a random
        lig_to_mutate = random.randint(0,2)
        child = octahedral_complex(self.ligands_dict)
        child.copy(self) # copies this parent
        n = len(self.ligands_dict)
        self.examine()
        print('I think this is 3x bidentate: ',self.three_bidentate,self.ax_dent)
        if (lig_to_mutate == 0):
            print("mutating equitorial")
            rand_ind = numpy.random.randint(low = 0,high = n)
            child.replace_equitorial([rand_ind])
        else:
            ready_flag = False
            while not ready_flag:
                new_ax_list = list()
                rand_ind = numpy.random.randint(low = 0,high = n)
                ax_ligand_properties  = self.ligands_dict[rand_ind][1]
                ax_dent = ax_ligand_properties[0]
                if (ax_dent == self.ax_dent):
                    if (lig_to_mutate == 1):
                        print("mutating axial 1 ")
                        new_ax_list = [rand_ind,self.ax_inds[1]]
                    elif (lig_to_mutate == 2):
                        print("mutating axial 2 ")
                        new_ax_list = [self.ax_inds[0],rand_ind]
                    child.ax_dent = 1
                    child.three_bidentate = False
                    ready_flag = True
                elif (ax_dent  == 2) and (self.ax_dent == 1):
                    ## here, we want to add a bidentate but 
                    ## need to swap the second ligand too
                    print("swapping both axial ")
                    new_ax_list = [rand_ind,rand_ind]
                    child.ax_dent = 2
                    child.three_bidentate = True
                    ready_flag =  True
                if set(new_ax_list) == set([0,0]):
                    ready_flag = False

            print("trying to add " + str(new_ax_list))
            child.replace_axial(new_ax_list)
        child._name_self()
        child.examine()
        return child
    def write_files(self):
        pass

    def choose_tolerances(self,size,f):
        if size <= 50:
            etol = '2.5e-6'
            grms = '2.5e-4'
            gmax = '5e-4'
            drms = '1.5e-3' 
            dmax = '3e-3'
        elif (size > 50) and (size <= 100):
            etol = '5e-6'
            grms = '5e-4'
            gmax = '1e-3'
            drms = '1.5e-3' 
            dmax = '5e-3'
        elif (size > 100) and (size <= 200):
            etol = '1e-5'
            grms = '1e-3'
            gmax = '5e-3'
            drms = '3e-3' 
            dmax = '6e-3'
        elif (size > 200):
            etol = '5e-5'
            grms = '1e-3'
            gmax = '5e-3'
            drms = '3e-3' 
            dmax = '6e-3'
        f.write('min_converge_e ' + etol +'\n')
        f.write('min_converge_grms ' + grms +'\n')
        f.write('min_converge_gmax ' + gmax +'\n')
        f.write('min_converge_drms ' + drms +'\n')
        f.write('min_converge_dmax ' + dmax +'\n')







    def _generate_geometery(self,prefix,spin):
        self._name_self()
        path_dictionary = setup_paths()
        rundirpath = get_run_dir()
        molsimpath = path_dictionary["molsimp_path"]
        ligloc_cont = True
        mol_name = prefix + self.name + "_" + str(spin)
        # defaults for ocathedrals
        geometry = "oct"
        cord =6
        self.check_self_loud()
        if (self.ax_dent ==1) and (self.geo == 'spy'):
            present_axlig = [i for i in self.ax_ligands if i!="vacant"];
            liglist = (str([str(element).strip("'[]'") for element in (self.eq_ligands)]  + [str(element).strip("'[]'") for element in present_axlig]).strip("[]")).replace("'","")
            ligloc = 1
            geometry = "spy"
            cord = 5
        elif self.ax_dent == 1:
           liglist = (str([str(element).strip("'[]'") for element in (self.eq_ligands)]  + [str(element).strip("'[]'") for element in self.ax_ligands]).strip("[]")).replace("'","")
           ligloc = 1
        elif self.ax_dent == 2:
           liglist = (str([str(element).strip("'[]'") for element in (self.eq_ligands)]  + [str(self.ax_ligands[0]).strip("'[]'")]).strip("[]")).replace("'", "") 
           ligloc = 0
        if self.oxidation_state == 2:
            ox_string = "II"
        elif self.oxidation_state == 3:
            ox_string = "III"
        ms_dump_path = path_dictionary["molsimplify_inps"] +  'ms_output.txt'

        jobpath = path_dictionary["job_path"]  + mol_name + '.in'
        ## check if already exists:
        geo_exists = os.path.isfile(path_dictionary["initial_geo_path"] + mol_name + '.xyz')
        in_exists = os.path.isfile(jobpath)
        if not (geo_exists and in_exists):
                print('generating '+mol_name,self.eq_ligands,self.ax_ligands)
                with open(ms_dump_path,'a') as ms_pipe:
        #            print('writing to output to pipe' + ms_dump_path)
                    call = [molsimpath + '/main.py','-core ' + self.core,'-lig ' +liglist,
                             '-rundir ' + rundirpath +'\n','-jobdir temp','-keepHs yes,yes,yes,yes,yes,yes',
                             '-coord '+str(cord),'-ligalign 1','-ligloc ' + str(ligloc),'-calccharge yes','-name ' +mol_name + '\n',
                             '-geometry ' + geometry,'-spin ' + str(spin),'-oxstate '+ ox_string,
                             '-qccode TeraChem','-runtyp minimize','-method UDFT','-mopac']
        #            print(call)
                    p2 = subprocess.call(call,stdout = ms_pipe)
                shutil.move(rundirpath + 'temp/' + mol_name + '.molinp', path_dictionary["molsimplify_inps"]+'/' + mol_name + '.molinp')
                shutil.move(rundirpath + 'temp/' + mol_name + '.xyz', path_dictionary["initial_geo_path"] +'/'+ mol_name + '.xyz')
                shutil.move(rundirpath + 'temp/' + mol_name + '.mop', path_dictionary["mopac_path"] +'/'+ mol_name + '.mop')

                size = harvest_size(path_dictionary['initial_geo_path']+ '/'+mol_name+'.xyz')
                with open(jobpath,'w') as newf:
                    newf.writelines("run minimize \n")
                    with open(rundirpath + 'temp/' + mol_name + '.in','r') as oldf: 
                        for line in oldf:
                            if not ("coordinates" in line) and (not "end" in line) and not ("scrdir" in line) and not("run" in line) and not ("maxit" in line):
                                newf.writelines(line)
                    newf.writelines("min_coordinates cartesian \n")
                    newf.writelines("scrdir scr/geo/" + mol_name + "\n")
                    self.choose_tolerances(size,newf) # fetch size-aware tolerances
                os.remove(rundirpath + 'temp/' + mol_name + '.in')
        return jobpath

#print("\n\n\n\n")
#susan = TM_complex('susan')
#susan.random_gen()
#susan.examine()
#susan.ge
#print("\n ~~~~~~~~~~~~ \n")
#jack = octahedral_complex(ligands_dict)
#jack.random_gen()
#jack.examine()
#print("\n ~~~~~~~~~~~~ \n")
#fred = susan.exchange_ligands(jack)
#fred.examine()
#print("\n ~~~~~~~~~~~~ \n")
#fred = susan.exchange_ligands(jack)
#fred.examine()
#print("\n ~~~~~~~~~~~~ \n")
#fred = fred.mutate()
#fred.examine()
#path_dictionary = setup_paths()
#susan.replace_equitorial([32])
#print("\n\n")
#susan.examine()
#susan.replace_axial([2,2])
#print("\n\n")
#susan.examine()

#susan.generate_geometery(prefix = "te_",spin = 3) 




