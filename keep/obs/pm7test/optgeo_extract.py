import sys
import os

def parse_xyz(f):

    ## written by Terry
    s = f.read()
    ss = s.splitlines()

    lines = [line.split() for line in ss]
    el = list()
    x = list()
    y = list()
    z = list()
    controlflag = 0
    for i in lines:
#        print(i)
        if controlflag:
            controlflag = 0
            el = list()
            x = list()
            y = list()
            z = list()
#           print('clearing')
        else:
            try:
                z.append(i[3])
                el.append(i[0])
                x.append(i[1])
                y.append(i[2])
            except:
                controlflag = 1
    natoms = len(el)
    crds = list()
    elements = list()
#    print(z)
    for i in range(natoms):
        crds.append([x[i],y[i],z[i]])
        elements.append(el[i])
    ret_dict =dict()
    ret_dict['coords'] = crds
    ret_dict['elements'] = elements
    return ret_dict  

## load parser module and explanatory text 
## argpase not on gib
#parser =argparse.ArgumentParser(description = 'Convert Terachem geo optim history into optimal  geo xyz')
#parser.add_argument("targetpath",help="path to the target file")
#parser.add_argument("newfilename",help="name of the new file")

#args = parser.parse_args()
########################################################################################################################
########################################################################################################################
## extract file name and extension
try:
    filename = os.path.basename(sys.argv[1])
    name, ext = os.path.splitext(filename)
except:
    print('Error: unable to strip file name. Please check that a valid path is given')
    sys.exit('Error: unable to strip file name. Please check that a valid path is given')
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#create  a packed argument:
## now call the correct funciton
print('Hello, running python code!')
print(sys.argv)
print("read from " + sys.argv[1])
print("targeting loop " + sys.argv[2])

with open(sys.argv[1],'r') as f:
        print('after read')
        ret_dict=parse_xyz(f)
print("before loop"+ sys.argv[2])

with open(sys.argv[2],'w') as f:
        print('In loop')
        f.write(str(len(ret_dict['elements']))+'\n')
        f.write('# extracted from Terachem optimization\n')
        for i,elements in enumerate(ret_dict['elements']):
                writebuffer = [elements] + ret_dict['coords'][i]
                f.write(' '.join(s for s in writebuffer) + '\n')
print("finished"+ sys.argv[2])



