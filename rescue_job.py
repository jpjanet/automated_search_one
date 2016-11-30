import shutil
import os
import sys
from prep_calc import *
print('hello, rescue scripted invoked!')
if (len(sys.argv) != 2):
    print('incorrect arguments to resuce script')
    print(sys.argv)
    raise ValueError()
## sanitize input
path = sys.argv[1].strip('\n')
print('path is  '+path)


## first, set up paths:

