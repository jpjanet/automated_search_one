import os
########################

def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
########################

def get_run_dir():
    rdir = "/home/jp/redox_search/" 
    return rdir
########################

def translate_job_name(job):
    #print('job  = ' + str(job))

    base = os.path.basename(job)
    base = base.strip("\n")
    #print('base  = ' + str(base))
    base_name = base.strip(".in")
    base_name = base_name.strip(".done")
    low_name = str(base_name).lower()
    ll = (str(base_name)).split("_")
    zoo = ll[0]
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
        print('spin assigned as ll[9]  = ' + str(spin) + ' on  ' +str(ll))
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
                       "molsimp_path"     : "/home/jp/redox_search/",
                       "mopac_path"     : working_dir + "mopac/"}

    for keys in path_dictionary.keys():
        ensure_dir(path_dictionary[keys])
    return path_dictionary
########################

def get_ligands():
        f_ligands_dict =     [['vacant',[1]], #0
                       ['water',[1]],#1
                       ['acetonitrile',[1]],#2
                       ['imidazole',[1]],#3
                       ['pph3',[1]],#4
                       ['pyr',[1]],#5
                       ['ammonia',[1]],#6
                       ['carbonyl',[1]],#7
                       ['thiane',[1]],#8
                       ['furan',[1]],#24
                       ['misc',[1]],#9
                       ['chloropyridine',[1]],#14
                       ['pisc',[1]],#10
                       ['bipy',[2]],#11
                       ['phen',[2]],#12
                       ['en',[2]],#13
                       ['dppe',[2]],#19
                       ['pme3',[1]],#25
                       ['phosphine',[1]],#26
                       ['tetrahydrofuran',[1]],#27
                       ['thiopyridine',[1]],#28
                       ['mebipyridine',[2]],#29
                       ['ethbipyridine',[2]],#30
                       ['ethOHbipyridine',[2]],#31
                       ['nitrobipyridine',[2]],#32
                       ['phosacidbipyridine',[2]],#33
                       ['sulfacidbipyridine',[2]],#34
                       ['phendione',[2]],#35
                       ['phenacac',[2]],#36
                       ['tbisc',[1]],#37
                       ['phenisc',[1]],#38
                       ['porphyrin',[4]],#15
                       ['phthalocyanine',[4]],#16
                       ['cyclam',[4]],#17
                       ['phenylcyc',[4]],#18
                       ['corrolazine',[4]],#20
                       ['cyclen',[4]],#21
                       ['tbutylcyclen',[4]],#22
                       ['cyanoaceticporphyrin',[4]]]

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
                       ['acac',[2]],#23
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
                       ['tetrahydrofuran',[1]],#50
                       ['thiopyridine',[1]],#51
                       ['mebipyridine',[2]],#52
                       ['ethbipyridine',[2]],#53
                       ['ethOHbipyridine',[2]],#54
                       ['nitrobipyridine',[2]],#55
                       ['phosacidbipyridine',[2]],#56
                       ['sulfacidbipyridine',[2]],#57
                       ['phendione',[2]],#58
                       ['benzenethiol',[1]],#59
                       ['benzenedithiol',[2]],#60
                       ['quinoxalinedithiol',[2]], #61
                       ['phenacac',[2]],#62
                       ['tbisc',[1]],#63
                       ['phenisc',[1]],#64
                       ] #25}
        
        return ligands_dict
########################

def get_metals():
        metals_list = ['cr','mn','fe','co']
        return metals_list
########################

def spin_dictionary():
    metal_spin_dictionary = {'co':{2:[2,4],3:[1,5]},
                              'cr':{2:[3,5],3:[2,4]},
                              'fe':{2:[1,5],3:[2,6]},
                              'mn':{2:[2,6],3:[1,5]}}
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

def set_outstanding_jobs(list_of_jobs):
    path_dictionary = setup_paths()
    path = path_dictionary['job_path']
    ensure_dir(path)
    print('settting jobs to be ' + str(list_of_jobs))
    with open(path + '/outstanding_job_list.txt', 'w') as f:
        for jobs in list_of_jobs:
            f.write(jobs.strip("\n") + "\n")
    print('written\n')
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
def find_submitted_jobs():
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
    return converged_job_dictionary
########################
def update_converged_job_dictionary(jobs,status):
        path_dictionary = setup_paths()
        converged_job_dictionary = find_converged_job_dictionary()
        converged_job_dictionary.update({jobs:status})
        if status == 2:
                print(' wrtiting job as s2 ' + str(jobs))
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

######################

