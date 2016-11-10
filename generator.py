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
                self.genes =  dict()
                self.gene_convergence_dictionary = dict()
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
                        this_unique_name = this_complex.name
                        print('this this_unique_name ', this_unique_name)
                        if not this_unique_name in self.unique_complex_dictionary.keys():
                            ## we can accept this complex
                            self.genes[ID] = this_unique_name
                            self.unique_complex_dictionary[ID] = this_complex
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
                print(this_complex.ax_inds)
                this_unique_name = this_complex.name
                print('this this_unique_name ', this_unique_name)
                if not this_unique_name in self.unique_complex_dictionary.keys():
                   ## we can accept this complex
                   self.genes[ID] = this_unique_name
                   self.unique_complex_dictionary[ID] = this_complex
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
            print('mox is ' + str(metal) + ' '+ str(this_ox))
            return this_ox[0]
        def write_state(self):
                ## first write genes to path
                state_path = self.base_path_dictionary["state_path"] +"current_complexes.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'w').close()
                else:   ## backup state data
                        shutil.copyfile(state_path,self.base_path_dictionary["state_path"] +"current_complexes.csv.bcp")
                emsg = write_dictionary(self.genes,state_path)
                ## second write live info to base directory
                state_path = self.base_path_dictionary["state_path"] +"current_status.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'w').close()
                emsg = write_dictionary(self.status_dictionary,state_path)
                if emsg:
                        print(emsg)

                ## third,  write gene-fitness info to path
                state_path = self.base_path_dictionary["state_path"] +"complex_convergence.csv"
                if not os.path.isfile(state_path):
                        open(state_path,'w').close()
                emsg = write_dictionary(self.gene_convergence_dictionary,state_path)
                if emsg:
                        print(emsg)


        def report(self):
            print('My name is ' + self.name)
            print('there are '+ str(self.total_counter)+ ' unique M-L combinations:')
            for keys in self.unique_complex_dictionary.keys():
                print(keys,self.unique_complex_dictionary[keys].name)
                self.name_assemble(keys,self.unique_complex_dictionary[keys].name)
        def name_assemble(self,ID,name):
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
                state_path = self.base_path_dictionary["state_path"] +"current_complexes.csv"
                emsg,complex_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                for ID in complex_dict.keys():
                    this_complex = self.name_assemble(int(ID),complex_dict[ID])
                    self.genes[int(ID)] = this_complex.name 
                    self.unique_complex_dictionary[int(ID)] = this_complex

                ## third,  read gene-fitness info to path
                state_path = self.base_path_dictionary["state_path"] +"complex_convergence.csv"
                emsg,fit_dict = read_dictionary(state_path)
                if emsg:
                        print(emsg)
                self.complex_convergence_dictionary = fit_dict

        def assess_convergence(self,completed_jobs):
            ## loop all over genes in the pool and the selected set
            self.outstanding_jobs = dict()
            for uniqkeys in self.genes.keys(): ## keys here are ID numbers
                this_names = self.genes[uniqkeys]
                ## see if this gene is in the fitness dictionary
                if this_names in completed_jobs:
                    pass
                    ## here, the job 
                else:
                    self.outstanding_jobs.update({uniqkeys:self.unique_complex_dictionary[uniqkeys]})
            logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) + " : " 
                       +  str(len(self.outstanding_jobs.keys())) + " calculations to be completed")
            print('length of jobskeys',len(self.outstanding_jobs.keys()))
            if (len(self.outstanding_jobs.keys()) ==0):
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) 
                               + " all jobs completed ")
            else:
                self.job_dispatcher()


        def job_dispatcher(self):
                jobpaths = list()
                msd = spin_dictionary()
                spin_lab = ['LS','HS']
                for keys in self.outstanding_jobs.keys():
                        #each key is a job name #
                        jobs = self.outstanding_jobs[keys]
                        this_metal = jobs.core
                        this_ox = jobs.oxidation_state
                        these_spins = msd[this_metal][this_ox]
                        print(these_spins)
                        for i in [0,1]:
                                job_prefix =  str(keys) + '_' 
                                ## generate HS/LS
                                this_jobname = jobs._generate_geometery(prefix = job_prefix, spin = these_spins[i])
                add_jobs(self.base_path_dictionary['job_path'],jobpaths)


        def select_best_genes(self):
                ## first write genes to path
                summary_path = self.current_path_dictionary["state_path"] +"all_genes.csv"
                outcome_list = list()
                npool  =  self.status_dictionary["npool"]
                mean_fitness = 0
                for keys in self.genes.keys():
                    outcome_list.append((keys,self.genes[keys],float(self.gene_fitness_dictionary[self.genes[keys]])))
                    logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                               ": Gen " + str(self.status_dictionary['gen']) 
                             + " fitness is = " +  str(float(self.gene_fitness_dictionary[self.genes[keys]])))
 
                outcome_list.sort(key=lambda tup: tup[2], reverse = True)

                full_size = len(outcome_list)

                if not os.path.isfile(summary_path):
                       open(summary_path,'a').close()
                emsg = write_summary_list(outcome_list,summary_path)
                self.genes = dict()
                self.gene_compound_dictionary = dict()
                for i in range(0,npool):
                        self.genes[i] = outcome_list[i][1]
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.encode(self.genes[i])
                        self.gene_compound_dictionary[i] = this_complex
                        mean_fitness += float(outcome_list[i][2])
                mean_fitness = mean_fitness/npool # average fitness
                self.status_dictionary.update({'mean_fitness':mean_fitness})
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                       ": Gen " + str(self.status_dictionary['gen']) 
                     + " complete, mean_fitness = " +  str(mean_fitness))
        def advance_generation(self):
                ## advance counter
                self.status_dictionary['gen'] +=1
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                       ": Gen " + str(self.status_dictionary['gen']-1) 
                     + " advancing to Gen " +  str(self.status_dictionary['gen']))
                self.status_dictionary['ready_to_advance'] = False
                self.current_path_dictionary = advance_paths(self.base_path_dictionary,self.status_dictionary['gen'])

                npool  =  self.status_dictionary["npool"]
                ncross =  self.status_dictionary["ncross"]
                pmut   =  self.status_dictionary["pmut"]

                ## generation selected set
                selected_genes = dict()
                selected_compound_dictionary = dict()
                number_selected = 0
                ## populate selected pool
                while number_selected < npool:
                        this_int = random.randint(0,npool -1)
                        this_barrier = random.uniform(0,1)
                        this_gene = self.genes[this_int]
                        if self.gene_fitness_dictionary[this_gene] > this_barrier:
                                selected_genes[number_selected + npool] = this_gene
                                number_selected += 1
                ## populate compound list
                for keys in selected_genes.keys():
                        genes = selected_genes[keys]
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.encode(genes)
                        selected_compound_dictionary[keys] = this_complex
                ## now perfrom ncross exchanges
                number_of_crosses = 0
                while number_of_crosses < ncross:
                        these_partners = random.sample(range(npool,(2*npool - 1)),2)
                        keep_axial = selected_compound_dictionary[these_partners[0]]
                        keep_equitorial = selected_compound_dictionary[these_partners[1]]
                        old_genes = [selected_genes[key] for key in these_partners]
                        new_complex_1 = keep_axial.exchange_ligands(keep_equitorial,True)
                        new_complex_2 = keep_equitorial.exchange_ligands(keep_axial,True)
                        new_gene_1 = new_complex_1.name
                        new_gene_2 = new_complex_2.name
                        selected_genes[these_partners[0]] = new_gene_1
                        selected_compound_dictionary[these_partners[0]] = new_complex_1
                        selected_genes[these_partners[1]] = new_gene_2
                        selected_compound_dictionary[these_partners[1]] = new_complex_2
                        new_genes = [selected_genes[key] for key in these_partners]

                        number_of_crosses +=1
                        logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                               ":  Gen " + str(self.status_dictionary['gen'])
                               + " crossing " + str(these_partners) + " " +
                              str(old_genes) + " -> " + str(new_genes)  )

                ## mutate
                for keys in selected_genes.keys():
                        does_mutate = random.uniform(0,1)
                        if does_mutate < pmut:
                                print("\n")
                                old_gene = selected_genes[keys]
                                mutant = selected_compound_dictionary[keys].mutate()
                                selected_compound_dictionary[keys] = mutant
                                selected_genes[keys] = selected_compound_dictionary[keys].name
                                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                                       ":  Gen " + str(self.status_dictionary['gen'])
                                       + " mutating " + str(keys) + ": "  + old_gene +  " -> " + mutant.name) 

                ## merge the lists 
                self.genes.update(selected_genes)
                self.gene_compound_dictionary.update(selected_compound_dictionary)



try_pop = complex_generator('Jenny')
try_pop.populate_random(10)
try_pop.write_state()
try_pop.report()

second_pop=complex_generator('Justin')
second_pop.read_state()
#for i in range(1,47):
#    print(i)
#    second_pop.populate_lig_combo([[i],[1, 1]])
#second_pop.populate_lig_combo([[41],[1, 1]])

second_pop.report()
second_pop.assess_convergence(dict())


