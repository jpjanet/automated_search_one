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
    ID,low_name,base_name,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat,gene = translate_job_name(job)

    base_name = os.path.basename(job).strip('in')
    if sub_num > 1:
        print(' start rescue')
        ## run rescue and analysis
#        rescue_cmd_str = './gibraltar_rescue_.sh ' + job
#        p_res = subprocess.Popen(cmd_str,shell=True,stdout=subprocess.PIPE)
    ## could call different script if resub? currently only calls the same
    cmd_str ='qsub -j y -N  ' +'r_'+ str(ID) + '_'+str(spin_cat) + ' ' +get_run_dir() + 'gibraltar_wrap_auto.sh ' + job
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
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    number_live_jobs = len(live_job_dictionary.keys())
    ## set of jobs to dispatch
    joblist  = get_outstanding_jobs()
    logger(path_dictionary['state_path'],str(datetime.datetime.now())
           + " number of calculations to be completed =   " + str(len(joblist)))

    sub_count = 0;
    resub_count = 0;
    lmax = 30 #number of live jobs
    if number_live_jobs < lmax:
        for jobs in joblist:
            jobs = jobs.strip("\n")
            if (not (jobs in live_job_dictionary.keys())) and (len(jobs.strip('\n')) != 0 ) and (number_live_jobs < lmax): ## check the job isn't live
                print(jobs,'is not live....')
#                print('is it in sub:' + str(submitted_job_dictionary.keys()))
                if not (jobs in submitted_job_dictionary.keys()):
                    ## launch
                    submitted_job_dictionary.update({jobs:1})
                    ## submit job to queue
                    job_id = launch_job(jobs,1)
                    sub_count += 1
                    number_live_jobs += 1
                    print('updating LJD with :',job_id,jobs)
                    live_job_dictionary.update({jobs:job_id})
                else: # job is a resubmission 
                    number_of_attempts = submitted_job_dictionary[jobs]
                    print('number of attempts = '+ str(number_of_attempts))
                    if (int(number_of_attempts) <= 12):
                        ## relaunch  
                        submitted_job_dictionary.update({jobs: (int(number_of_attempts)+1)})
                        job_id = launch_job(jobs,int(number_of_attempts) + 1)
                        number_live_jobs += 1
                        resub_count += 1
                        print('(resub: '+str(resub_count)+ ' )updating LJD with :' + str(job_id) + ' ' + str(jobs))
                        live_job_dictionary.update({jobs:job_id})

                    else: # give up on this job 
                        logger(path_dictionary['state_path'],str(datetime.datetime.now())
                           + " Giving up on job : " + str(jobs) + ' with '+ str(number_of_attempts) + ' attempts')
                        update_converged_job_dictionary(jobs,2) # mark job as abandoned 
            else:
                print('job is live or empty or queue is full')
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
    submitted_job_dictionary = find_submitted_jobs()
    ## live jobs:
    live_job_dictionary = find_live_jobs()
    ## check the status of current jobs
    joblist  = get_outstanding_jobs()
    all_runs = dict()
    jobs_complete = 0
    for jobs in joblist:
        if (jobs not in live_job_dictionary.keys()) and ((len(jobs.strip('\n'))!=0)) and (jobs not in submitted_job_dictionary.keys()):
 #           print(jobs)
            this_run = test_terachem_go_convergence(jobs)
            print("is this run conv: " + str(this_run.converged)+ ' at ' + str(jobs))
            if this_run.converged == True:
                jobs_complete += 1
                update_converged_job_dictionary(jobs,1) # record converged 
                remove_from_outstanding_jobs(jobs) # take out of pool
        else:
                print('check_all thinks we are live: '+ str(jobs))
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
 #           print('job_id and status : ',this_job_id,this_status)
            if this_status:
                counter += 1
#                print('recording as live:'+ str(jobs) +str(this_job_id))
                live_job_dictionary.update({jobs:this_job_id})
            else:
                    print('del job: '+ str(jobs))
                    del live_job_dictionary[jobs]
    write_dictionary(live_job_dictionary,
                     path_dictionary["job_path"]+"/live_jobs.csv")
    return counter
########################
def remove_from_outstanding_jobs(job):
        current_outstanding = get_outstanding_jobs()
 #       print('job ' + str(job)+ ' is  complete, removing..')
 #       print(current_outstanding)
        if job in current_outstanding:
                current_outstanding.remove(job)
                print('\n job found and removed')
#        print('\n****************************\n')
#        print(current_outstanding)
  #      sad
        set_outstanding_jobs(current_outstanding)
########################




