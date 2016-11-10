import os

def ensure_dir(dir_path):
    if not os.path.exists(dir_path):
        print('creating' + dir_path)
        os.makedirs(dir_path)
def get_run_dir():
    rdir = "/home/jp/Dropbox/Main/optimal_mol_des/automated_search_one/" 
    return rdir

def translate_job_name(job):
    base = os.path.basename(job)
    base = base.strip("\n")
    base_name = base.strip(".in")
    ll = (str(base)).split("_")
    slot = ll[4]
    gen = int(ll[1])
    gene = str(ll[4]+"_"+ll[5]+"_"+ll[6])
    spin = int(ll[7].rstrip(".in"))
    return gen,slot,gene,spin,base_name


def setup_paths():
    working_dir = get_run_dir()
    print('working dir ' + working_dir)
    path_dictionary = {"out_path"     : working_dir + "outfiles/",
                   "job_path"         : working_dir + "jobs/",
                   "done_path"        : working_dir + "completejobs/",
                   "initial_geo_path" : working_dir + "initial_geo/",
                   "optimial_geo_path": working_dir + "optimized_geo/",
                   "state_path"       : working_dir + "statespace/",
                   "molsimplify_inps" : working_dir + "ms_inps/",
                   "infiles"          : working_dir + "infiles/",
                   "molsimp_path"     : working_dir + "molSimplify/",
                   "mopac_path"     : working_dir + "mopac/"}

    for keys in path_dictionary.keys():
        ensure_dir(path_dictionary[keys])
    return path_dictionary
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
                       ['acac',[2]],#23
                       ['en',[2]],#24
                       ['tbuc',[2]],#25
#                       ['bpabipy',[4]],#26
                       ['porphyrin',[4]],#27
                       ['phthalocyanine',[4]],#28
                       ['cyclam',[4]],#29
                       ['phenylcyc',[4]],
                       ['dppe',[2]],#30
                       ['corrolazine',[4]],#31
                       ['iminodiacetic',[2]],#32
                       ['diaminomethyl',[2]],#33
                       ['cyclen',[4]],#34
                       ['tbutylcyclen',[4]],#35
                       ['cyanoaceticporphyrin',[4]],#36
                       ['dicyanamide',[1]],#37
#                       ['tcnoet',[1]],#38
#                       ['tcnoetOH',[1]],#39
                       ['methanethiol',[1]],#40
                       ['ethanethiol',[1]],#41
                       ['tbutylthiol',[1]],#42
                       ['pme3',[1]],#43
                       ['tricyanomethyl',[1]],#44
                       ['pentacyanopentadienide',[1]],#45
                       ['phosphine',[1]],#46
                       ['propdiol',[2]]
                       ] #25}
        
    return ligands_dict
def get_metals():
        metals_list = ['Cr','Mn','Fe','Co','Ni']
        return metals_list
def spin_dictionary():
    metal_spin_dictionary = {'Co':{2:[2,4],3:[1,5]},
                              'Cr':{2:[3,5],3:[2,4]},
                              'Fe':{2:[1,5],3:[2,6]},
                              'Mn':{2:[2,6],3:[3,5]},
                              'Ni':{2:[1,3]}}
    return metal_spin_dictionary

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
def logger(path, message):
    ensure_dir(path)
    with open(path + '/log.txt', 'a') as f:
        f.write(message + "\n")
def add_jobs(path,list_of_jobs):

    ensure_dir(path)
    all_jobs = load_jobs(path)
    for items in list_of_jobs:
            if not item in all_jobs:
                    all_jobs.append(item)
    with open(path + '/current_jobs.txt', 'w') as f:
        for jobs in all_jobs:
            f.write(jobs + "\n")
def load_jobs(path):
    ensure_dir(path)
    list_of_jobs = list()
    if os.path.exists(path + '/current_jobs.txt'):
        with open(path + '/current_jobs.txt', 'r') as f:
            for lines in f:
                list_of_jobs.append(lines)
    return list_of_jobs
def harvest_size(filename):
    size = 0;
    if not os.path.exists(filename):
        print('error, file not found')
    else:
        with open(filename,'r') as f:
            ff = f.readlines()
            print('size is ' + str(ff[0]))
            size = int(ff[0])
    return size


