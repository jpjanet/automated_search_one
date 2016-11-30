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
                        this_gene = self.gene_id_dictionary[this_int]
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
                self.gene_id_dictionary.update(selected_genes)
                self.gene_compound_dictionary.update(selected_compound_dictionary)


        def select_best_genes(self):
                ## first write genes to path
                summary_path = self.current_path_dictionary["state_path"] +"all_genes.csv"
                outcome_list = list()
                npool  =  self.status_dictionary["npool"]
                mean_fitness = 0
                for keys in self.genes.keys():
                    outcome_list.append((keys,self.genes[keys],float(self.gene_fitness_dictionary[self.gene_id_dictionary[keys]])))
                    logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                               ": Gen " + str(self.status_dictionary['gen']) 
                             + " fitness is = " +  str(float(self.gene_fitness_dictionary[self.gene_id_dictionary[keys]])))
                outcome_list.sort(key=lambda tup: tup[2], reverse = True)

                full_size = len(outcome_list)

                if not os.path.isfile(summary_path):
                       open(summary_path,'a').close()
                emsg = write_summary_list(outcome_list,summary_path)
                self.gene_id_dictionary = dict()
                self.gen_complex_dictionary = dict()
                for i in range(0,npool):
                        self.gene_id_dictionary[i] = outcome_list[i][1]
                        this_complex = octahedral_complex(self.ligands_list)
                        this_complex.encode(self.gene_id_dictionary[i])
                        self.unique_complex_dictionary[i] = this_complex
                        mean_fitness += float(outcome_list[i][2])
                mean_fitness = mean_fitness/npool # average fitness
                self.status_dictionary.update({'mean_fitness':mean_fitness})
                logger(self.base_path_dictionary['state_path'],str(datetime.datetime.now()) +
                       ": Gen " + str(self.status_dictionary['gen']) 
                     + " complete, mean_fitness = " +  str(mean_fitness))


