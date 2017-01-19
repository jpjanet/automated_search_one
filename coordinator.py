import glob
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
from generator import *

def initialize_complexes(n):
        path_dictionary = setup_paths()
        population = complex_generator('current') 
        population.populate_random(n)
        population.write_state()
        logger(population.base_path_dictionary['state_path'],str(datetime.datetime.now())
               + ": <new pop>  with : " + str(population.status_dictionary['total_counter']) + ' unique species ')

        return population

def wake_up_routine():
        ## set up environment:        
        path_dictionary = setup_paths()
        ## initialize class
        new_pop = complex_generator('current')
        ## read in info
        new_pop.read_state()
        logger(new_pop.base_path_dictionary['state_path'],str(datetime.datetime.now())
               + ": <resuming>  ")
        ## assess current fitness
        new_pop.assess_convergence()
        state_path = new_pop.base_path_dictionary["state_path"] +"converged_job_dictionary.csv"
        converged_job_dictionary = read_dictionary(state_path)
        outstanding_jobs = get_outstanding_jobs()
        outstanding_count  = len(outstanding_jobs)
        print('outstanding_count: '+ str(outstanding_count))
        if outstanding_count < 20:
            print('adding jobs')
            new_pop.populate_random(math.ceil(10-outstanding_count)/2)
        new_pop.write_state()
        logger(new_pop.base_path_dictionary['state_path'],str(datetime.datetime.now())
               + ": number of jobs outstanding is  "+str(outstanding_count))
        return(new_pop)


#pop = initialize_complexes(2)
#pop.report()
print('******************************')
#new_pop = wake_up_routine()
#new_pop.report()
#new_pop = wake_up_routine()
#new_pop.report()
