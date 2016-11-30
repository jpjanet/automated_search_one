import glob
import operator
import datetime
import math
import numpy
import subprocess
import argparse
import os
import random
import shutil
from prep_calc import *
from complex_classes import *

class complex_generator:
        def __init__(self,name):
                path_dictionary = setup_paths()
                ligands_list = get_ligands()
                self.base_path_dictionary = path_dictionary
                self.name = name
                self.gene_id_dictionary =  dict()
                self.complex_convergence_dictionary = dict()
                self.ligands_list = ligands_list
                self.status_dictionary = dict()
                self.unique_complex_dictionary = dict()
                self.total_counter = 0;
                self.configure_status(0)
        def configure_status(self,total_counter):
            self.status_dictionary.update({'total_counter':total_counter})
            self.total_counter = total_counter
        def populate_random(self,n):
                ## clear the pool
                ### fill the pool with random structures 
                this_counter  = 0
                while this_counter < n:
                        this_metal = self.random_metal()
                        ID = (1000+self.total_counter)
                        this_ox = self.random_ox(this_metal)
                        this_complex = TM_complex(ID,this_metal,this_ox)
                        this_complex.random_gen()
                        this_gene = this_complex.name
                        print('this this_unique_name ', this_gene)
                        if not this_gene in self.unique_complex_dictionary.keys():
                            ## we can accept this complex
                            self.gene_id_dictionary[ID] = this_gene
                            self.unique_complex_dictionary[this_gene] = this_complex
                            this_counter = this_counter + 1
                            self.total_counter = self.total_counter + 1 
                self.configure_status(self.total_counter)
        def populate_lig_combo(self,ligs):
                ## clear the pool
                this_metal = self.random_metal()
                ID = (1000+self.total_counter)
                this_ox = self.random_ox(this_metal)
                this_complex = TM_complex(ID,this_metal,this_ox)
                this_complex.random_gen()
                this_complex.replace_equitorial(ligs[0])
                this_complex.replace_axial(ligs[1])
#                print(this_complex.ax_inds)
                this_gene = this_complex.name
                print('this this_unique_name ', this_gene)
                if not this_gene in self.unique_complex_dictionary.keys():
                   ## we can accept this complex
                   self.gene_id_dictionary[ID] = this_gene
                   self.unique_complex_dictionary[this_gene] = this_complex
                   self.total_counter = self.total_counter + 1 
                self.configure_status(self.total_counter)
                print('\n')
        def populate_metal_ox_lig_combo(self,metal,ox,ligs):
                ## fetch metal and oxidation state
                this_metal= metal
                this_ox = ox
                print('metal is ' + str(this_metal))
                print('this ox is ' + str(this_ox))
                # assign ID + generate
                ID = (1000+self.total_counter)
                this_complex = TM_complex(ID,this_metal,this_ox)
                this_complex.random_gen()
                this_complex.replace_equitorial(ligs[0])
                this_complex.replace_axial(ligs[1])
                this_gene = this_complex.name
                print('this this_unique_name ', this_gene)
                if not this_gene in self.unique_complex_dictionary.keys():
                   ## we can accept this complex
                   self.gene_id_dictionary[ID] = this_gene
                   self.unique_complex_dictionary[this_gene] = this_complex
                   self.total_counter = self.total_counter + 1
                self.configure_status(self.total_counter)
                print('\n')

        def random_metal(self):
            metals_list = get_metals()
            metal= random.sample(metals_list,1)
            return metal[0]
        def random_ox(sel,metal):
            allowed_states_dictionary = spin_dictionary()
            aha =  allowed_states_dictionary[metal]
            this_ox = random.sample(aha.keys(),1)
#            print('mox is ' + str(metal) + ' '+ str(this_ox))
            return this_ox[0]
        def write_state(self):
                ## first write genes to path
                state_path = self.base_path_dictionary["state_path"] + "gene_id_dictionary.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'w').close()
                else:   ## backup state data
                        shutil.copyfile(state_path,self.base_path_dictionary["state_path"] +"gene_id_dictionary.csv.bcp")
                emsg = write_dictionary(self.gene_id_dictionary,state_path)
                ## second write live info to base directory
                state_path = self.base_path_dictionary["state_path"] +"current_status.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'w').close()
                emsg = write_dictionary(self.status_dictionary,state_path)
                if emsg:
                        print(emsg)

                ## third,  write gene-fitness info to path
                state_path = self.base_path_dictionary["state_path"] +"complex_convergence_dictionary.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'w').close()
                emsg = write_dictionary(self.complex_convergence_dictionary,state_path)
                if emsg:
                        print(emsg)


        def report(self):
            print('My name is ' + self.name)
            print('there are '+ str(self.total_counter)+ ' unique M-L combinations:')
            for keys in self.unique_complex_dictionary.keys():
#                print(keys,self.unique_complex_dictionary[keys].name)
                self.name_assemble(keys,self.unique_complex_dictionary[keys].name)
        def name_assemble(self,ID,name):
            ## convert gene to complex
            name_split = name.split('_')
            metal = name_split[0]
            ox = int(name_split[1])
            eqlist =[int(name_split[3])]
            axlist =[int(name_split[5]),int(name_split[7])]
            this_complex = TM_complex(ID,metal,ox)
            this_complex.encode([eqlist[0],axlist[0],axlist[1]])
            return this_complex

        def read_state(self):
                ## first read gve info from base directory
                state_path = self.base_path_dictionary["state_path"] +"current_status.csv"
                emsg,read_dict = read_dictionary(state_path)
                print(read_dict)
                if emsg:
                        print(emsg)
                self.configure_status(int(read_dict["total_counter"]))
                ## next read  genes from path
                state_path = self.base_path_dictionary["state_path"] +"gene_id_dictionary.csv"
                emsg,complex_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                for ID in complex_dict.keys():
                    this_complex = self.name_assemble(int(ID),complex_dict[ID])
                    self.gene_id_dictionary[int(ID)] = this_complex.name 
                    self.unique_complex_dictionary[this_complex.name] = this_complex

                ## third,  read gene-fitness info to path
                state_path = self.base_path_dictionary["state_path"] +"complex_convergence_dictionary.csv"
                emsg,fit_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                self.complex_convergence_dictionary = fit_dict

        def assess_convergence(self):
            ## complex convergence dictionary is read in with self 
            ## get job_convergence_dictionary:
            state_path = self.base_path_dictionary["state_path"] +"converged_job_dictionary.csv"
            converged_job_dictionary = read_dictionary(state_path)
            ## loop all over all comlpexes 
            self.ID_for_dispatch = list()
            print('unique_complex_dictionary keys:',self.unique_complex_dictionary.keys())
            for uniqkeys in self.gene_id_dictionary.keys(): ## keys here are ID numbers
                already_converged = False # test if need to submit jobs for this gene
                this_name = self.gene_id_dictionary[uniqkeys] ## these are the genes
                ## see if this gene is in the fitness dictionary
                if this_name in self.complex_convergence_dictionary.keys():
                        if int(self.complex_convergence_dictionary[this_name]) != 0:
                                already_converged = True
                if not already_converged:
                        self.ID_for_dispatch.append(uniqkeys)
            logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + " : " 
                       +  str(len(self.ID_for_dispatch)) + " calculations to be completed")
            print('length of dispatch list: ' + str(2*len(self.ID_for_dispatch)))
            if (len(self.ID_for_dispatch) ==0):
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) 
                               + " all jobs completed ")
            else:
                self.job_dispatcher()
        def job_dispatcher(self):
                jobpaths = list()
                msd = spin_dictionary()
                for ID in self.ID_for_dispatch:
                        #each key is a job name #
                        this_name = self.gene_id_dictionary[ID]
                        this_complex = self.unique_complex_dictionary[this_name]
                        this_metal = this_complex.core
                        this_ox = this_complex.oxidation_state
                        these_spins = msd[this_metal][this_ox]
#                        print(these_spins)
                        for i in [0,1]:
                                job_prefix =  str(ID) + '_' 
                                ## generate HS/LS
                                this_jobname = this_complex._generate_geometery(prefix = job_prefix, spin = these_spins[i])
                                jobpaths.append(this_jobname)
                set_outstanding_jobs(self.base_path_dictionary['job_path'],jobpaths)



try_pop = complex_generator('Jenny')
#try_pop.populate_random(4)
try_pop.write_state()
try_pop.report()

#second_pop=complex_generator('Justin')
#second_pop.read_state()

metals_list = get_metals()
for metal_ind in [0,1,2,3,4] :
    this_metal= metals_list[metal_ind]
    allowed_states_dictionary = spin_dictionary()
    allowed_ox =  allowed_states_dictionary[this_metal]
    for ox in allowed_ox:
        try_pop.populate_metal_ox_lig_combo(this_metal,ox,[[2],[2,2]])
#    second_pop.populate_lig_combo([[i],[1, 1]])
try_pop.assess_convergence()
#second_pop.populate_lig_combo([[27],[1, 1]])
#second_pop.populate_lig_combo([[37],[1, 1]])
#second_pop.populate_lig_combo([[35],[1, 1]])
#second_pop.populate_lig_combo([[36],[1, 1]])
#print(second_pop.status_dictionary['total_counter'])
#second_pop.populate_lig_combo([[30],[1, 1]])
#second_pop.write_state()
#print(second_pop.status_dictionary['total_counter'])

#second_pop.populate_lig_combo([[1],[39,39]])
#second_pop.populate_lig_combo([[45],[45, 45]])
#second_pop.populate_lig_combo([[45],[1, 1]])



#second_pop.report()
#second_pop.assess_convergence(dict())


