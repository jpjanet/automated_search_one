import os
########################

def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        print('creating' + dir_path)
        os.makedirs(dir_path)
########################

def get_run_dir():
    rdir = "/home/jp/Dropbox/Main/optimal_mol_des/automated_search_one/" 
    return rdir
########################

def translate_job_name(job):
    base = os.path.basename(job)
    base = base.strip("\n")
    base_name = base.strip(".in")
    base_name = base_name.strip(".done")
    print(base_name)
    low_name = str(base_name).lower()
    ll = (str(base_name)).split("_")
    ID = int(ll[0])
    metal = ll[1].lower()
    ox = int(ll[2])
    eqlig = int(ll[4])

    axlig1 = int(ll[6])
    axlig2 = int(ll[8])
    spin = int(ll[9])
    metal_list = get_metals()
    metal_ind = metal_list.index(metal)
    metal_spin_dictionary  = spin_dictionary()
    these_states = metal_spin_dictionary[metal][ox]
    if spin == these_states[0]:
            spin_cat = 'LS'
    elif spin == these_states[1]:
            spin_cat = 'HS'
    else:
        print('critical erorr, unknown spin: '+ str(spin))
    gene = str(metal_ind) + '_' + str(ox) +'_'+ str(eqlig) + '_' + str(axlig1) + '_' + str(axlig2) 
    return ID,low_name,base_name,metal,ox,eqlig,axlig1,axlig2,spin,spin_cat,gene

########################

def setup_paths():
    working_dir = get_run_dir()
#    print('working dir ' + working_dir)
    path_dictionary = {"geo_out_path"     : working_dir + "geo_outfiles/",
                       "vertical_out_path"     : working_dir + "vertical_outfiles/",
                       "job_path"         : working_dir + "jobs/",
                       "done_path"        : working_dir + "completejobs/",
                       "initial_geo_path" : working_dir + "initial_geo/",
                       "progress_geo_path": working_dir + "prog_geo/",
                       "optimial_geo_path": working_dir + "optimized_geo/",
                       "state_path"       : working_dir + "statespace/",
                       "molsimplify_inps" : working_dir + "ms_inps/",
                       "infiles"          : working_dir + "infiles/",
                       "molsimp_path"     : "/home/jp/Dropbox/MyGit/molSimplify_dev/molSimplify/",
                       "mopac_path"     : working_dir + "mopac/"}

    for keys in path_dictionary.keys():
        ensure_dir(path_dictionary[keys])
    return path_dictionary
########################

def get_ligands():
    ligands_dict =     [['vacant',[1]], #0
                       ['thiocyanate',[1]], #1
                       ['chloride',[1]],#2
                       ['water',[1]],#3
                       ['acetonitrile',[1]],#4
                       ['ethyl',[1]],#5
                       ['imidazole',[1]],#6
                       ['nitro',[1]],#7
                       ['pph3',[1]],#8
                       ['pyr',[1]],#9
                       ['trifluoromethyl',[1]],#10
                       ['methanal',[1]],#11
                       ['benzene',[1]],#12
                       ['isothiocyanate',[1]],#13
                       ['ammonia',[1]],#14
                       ['cyanide',[1]],#15
                       ['carbonyl',[1]],#16
                       ['thiane',[1]],#17
                      ['misc',[1]],#18
                       ['pisc',[1]],#19
                       ['bipy',[2]],#20
                       ['phen',[2]],#21
                       ['ox',[2]],#22
                       ['acac',[2]],#h
                       ['en',[2]],#24
                       ['tbuc',[2]],#25
                       ['chloropyridine',[1]],#26
                       ['porphyrin',[4]],#27
                       ['phthalocyanine',[4]],#28
                       ['cyclam',[4]],#29
                       ['phenylcyc',[4]],#30
                       ['dppe',[2]],#31
                       ['corrolazine',[4]],#32
                       ['iminodiacetic',[2]],#33
                       ['diaminomethyl',[2]],#34
                       ['cyclen',[4]],#35
                       ['tbutylcyclen',[4]],#36
                       ['cyanoaceticporphyrin',[4]],#37
                       ['dicyanamide',[1]],#38
                       ['furan',[1]],#39
#                       ['tcnoetOH',[1]],#39
                       ['methanethiol',[1]],#40
                       ['ethanethiol',[1]],#41
                       ['tbutylthiol',[1]],#42
                       ['pme3',[1]],#43
                       ['tricyanomethyl',[1]],#44
                       ['pentacyanopentadienide',[1]],#45
                       ['phosphine',[1]],#46
                       ['propdiol',[2]],#47
                       ['amine',[1]],#48
                       ['uthiol',[1]],#49
                       ['tetrahydrofuran',[1]],
                       ['thiopyridine',[1]],
                       ['mebipyridine',[2]],
                       ['ethbipyridine',[2]],
                       ['ethOHbipyridine',[2]],
                       ['nitrobipyridine',[2]],
                       ['phosacidbipyridine',[2]],
                       ['sulfacidbipyridine',[2]],
                       ['phendione',[2]],
                       ['benzenethiol',[1]],
                       ['benzenedithiol',[2]],
                       ['quinoxalinedithiol',[2]],
                       ['phenacac',[2]],
                       ['tbisc',[1]],
                       ['phenisc',[1]],
                       ] #25}
        
    return ligands_dict
########################

def get_metals():
        metals_list = ['cr','mn','fe','co','ni']
        return metals_list
########################

def spin_dictionary():
    metal_spin_dictionary = {'co':{2:[2,4],3:[1,5]},
                              'cr':{2:[3,5],3:[2,4]},
                              'fe':{2:[1,5],3:[2,6]},
                              'mn':{2:[2,6],3:[3,5]},
                              'ni':{2:[1,3]}}
    return metal_spin_dictionary
########################

def write_dictionary(dictionary,path,force_append = False):
    emsg =  False
    if force_append:
        write_control = 'a'
    else:
       write_control = 'w' 
    try:
       with open(path,write_control) as f:
            for keys in dictionary.keys():
                f.write(str(keys).strip("\n") + ',' + str(dictionary[keys]) + '\n')
    except:
        emsg = "Error, could not write state space: " + path
    return emsg
########################

def write_summary_list(outcome_list,path):
    emsg =  False
    try:
        with open(path,'w') as f:
            for tups  in outcome_list:
                for items in tups:
                    f.write(str(items) + ',')
                f.write('\n')
    except:
        emsg = "Error, could not write state space: " + path
    return emsg
########################

def read_dictionary(path):
    emsg =  False
    dictionary = dict()
    try:
        with open(path,'r') as f:
            for lines in f:
                ll = lines.split(",")
                key = ll[0]
                value = ll[1].rstrip("\n")
                dictionary[key] = value
    except:
        emsg = "Error, could not read state space: " + path
    return emsg,dictionary
########################

def logger(path, message):
    ensure_dir(path)
    with open(path + '/log.txt', 'a') as f:
        f.write(message + "\n")
########################

def set_outstanding_jobs(path,list_of_jobs):
    ensure_dir(path)
    with open(path + '/outstanding_job_list.txt', 'w') as f:
        for jobs in list_of_jobs:
            f.write(jobs + "\n")
########################

def get_outstanding_jobs():
    path_dictionary = setup_paths()
    path = path_dictionary['job_path']
    ensure_dir(path)
    print('looking in ' +path + '/outstanding_job_list.txt')
    list_of_jobs = list()
    if os.path.exists(path + '/outstanding_job_list.txt'):
        print('found job list')
        with open(path + '/outstanding_job_list.txt', 'r') as f:
            for lines in f:
                list_of_jobs.append(lines)
    return list_of_jobs
########################
def find_submited_jobs():
    path_dictionary = setup_paths()
    if os.path.exists(path_dictionary["job_path"]+"/submitted_jobs.csv"):
        emsg,submitted_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/submitted_jobs.csv")
    else:
        submitted_job_dictionary = dict()

    return submitted_job_dictionary
########################
def find_live_jobs():
    path_dictionary = setup_paths()
    live_job_dictionary = dict()
    if os.path.exists(path_dictionary["job_path"]+"/live_jobs.csv"):
        emsg,live_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/live_jobs.csv")
    else:
       live_job_dictionary = dict()
    return live_job_dictionary
########################
def find_converged_job_dictionary():
    path_dictionary = setup_paths()
    converged_job_dictionary = dict()
    if os.path.exists(path_dictionary["job_path"]+"/converged_job_dictionary.csv"):
            emsg,converged_job_dictionary = read_dictionary(path_dictionary["job_path"]+"/converged_job_dictionary.csv")
    else:
       converged_job_dictionary = dict()
########################
def update_converged_job_dictionary(jobs,status):
        path_dictionary = setup_paths()
        converged_job_dictionary = find_converged_job_dictionary()
        converged_job_dictionary.update({jobs:status})
        write_dictionary(converged_job_dictionary,path_dictionary["job_path"]+"/converged_job_dictionary.csv")

########################

def harvest_size(filename):
    size = 0;
    if not os.path.exists(filename):
        print('error, file not found')
    else:
        with open(filename,'r') as f:
            ff = f.readlines()
        #    print('size is ' + str(ff[0]))
            size = int(ff[0])
    return size
########################

def translate_job_name(job_name):
        ll = str(job_name).strip('"''"')
        ll = ll.strip('\n')
        ID = int(ll[0])
        metal = ll[1]
        ox = int(ll[2])
        eq_lig = int(ll[4])
        ax1_lig = int(ll[6])
        ax2_lig = int(ll[8])
        spin = int(ll[9])

        return ID,metal,ox,eq_lig,ax1_lig,ax2_lig,spin
########################

