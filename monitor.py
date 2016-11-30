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
from coordinator import *
from process_scf import *
#######################
def launch_job(job,sub_num):
    ## code to submit to queue
    print('lauching ' + job + ' sub number: '+ str(sub_num))
    basename = os.path.basename(job)
    if sub_num > 1:
        print(' start rescue')
        ## run rescue and analysis
#        rescue_cmd_str = './gibraltar_rescue_.sh ' + job
#        p_res = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ## could call different script if resub? currently only calls the same
    cmd_str ='qsub -j y -N ' + str(base_name) + ' ' +get_run_dir() + 'gibraltar_wrap_auto.sh ' + job
    p_sub = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ll = p_sub.communicate()[0]
    ll =  ll.split()
    job_id = ll[2]
    return job_id
########################
def is_job_live(job_id):
    cmd_str = ('qstat -j '+ str(job_id))
    p1 = subprocess.Popen(cmd_str,shell=True,
                          stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    rt,ll = p1.communicate()
    verdict = True
    ll=ll.split('\n')
    for lines in ll:
            if str(lines).find('Following jobs do not exist:') != -1:
                    print('job ' + str(job_id) + ' is not live')
                    verdict = False
    return verdict
########################
def submit_outstanding_jobs():
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submited_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    ## set of jobs to dispatch
    joblist  = get_outstanding_jobs()
    sub_count = 0;
    resub_count = 0;
    for jobs in joblist:
        jobs = jobs.strip("\n")
        print('job is ' +str(jobs))
        if (not (jobs in live_job_dictionary.keys())) and (len(jobs.strip('\n')) != 0 ): ## check the job isn't live
            print(jobs,'is not live....')
            if not (jobs in submitted_job_dictionary.keys()):
                ## launch
                submitted_job_dictionary.update({jobs:1})
                ## submit job to queue
                job_id = launch_job(jobs,1)
                sub_count += 1
                print('updating LJD with :',job_id,jobs)
                live_job_dictionary.update({jobs:job_id})
            else: # job is a resubmission 
                number_of_attempts = submitted_job_dictionary[jobs]
                if (number_of_attempts <= 3):
                    ## relaunch  
                    submitted_job_dictionary.update({jobs: (number_of_attempts+1)})
                    job_id = launch_job(jobs,number_of_attempts + 1)
                    resub_count += 1
                    print('(resub: '+str(resub_count)+ ' )updating LJD with :' + str(job_id) + ' ' + str(jobs))
                    live_job_dictionary.update({jobs:job_id})

                else: # give up on this job 
                    logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " Giving up on job after 3 attempts: : " + str(jobs))
                    update_converged_job_dictionary(jobs,2) # mark job as abandoned 
        else:
            print('job is live or empty')
    write_dictionary(submitted_job_dictionary, path_dictionary["job_path"] + "/submitted_jobs.csv")
    write_dictionary(live_job_dictionary, path_dictionary["job_path"] + "/live_jobs.csv")
    logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " submitted  " + str(sub_count) +' new jobs and ' + str(resub_count) + ' resubs ')

    return joblist
########################
def check_all_current_convergence():
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    submitted_job_dictionary = find_submmited_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    ## check the status of current jobs
    joblist  = get_outstanding_jobs()
    print(joblist)
    all_runs = dict()
    print(live_job_dictionary.keys())
    jobs_complete = 0
    for jobs in joblist:
        if jobs not in live_job_dictionary.keys() and (len(jobs.strip('\n'))!=0):
            print(jobs)
            this_run = test_terachem_sp_convergence(jobs)
            print("is this run conv: " + str(this_run.converged))
            if this_run.converged == True:
                jobs_complete += 1
                update_converged_job_dictionary({jobs:1}) # record converged 
    return jobs_complete
########################
def check_queue_for_live_jobs():
    ## This function reads the live_jobs csv 
    ## and tests if the jobs are still live 
    ## set up environment:        
    path_dictionary = setup_paths()
    ## previously dispatched jobs:
    ## live jobs in on record:
    live_job_dictionary = find_live_jobs()
    ## set of jobs requested by the algorithm
    counter = 0
    for jobs in live_job_dictionary.keys():
            this_job_id = live_job_dictionary[jobs]
            this_status = is_job_live(this_job_id)
            print('job_id and status : ',this_job_id,this_status)
            if this_status:
                counter += 1
                print('recording as live:'+ str(jobs) +str(this_job_id))
                live_job_dictionary.update({jobs:this_job_id})
            else:
                    if jobs in live_job_dictionary.keys():
                            del live_job_dictionary[jobs]
    write_dictionary(live_job_dictionary,
                     path_dictionary["job_path"]+"/live_jobs.csv")
    return counter




