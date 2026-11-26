import os, sys, glob
import numpy as np
from subprocess import call

#====set code path
CODE_PATH = '/Users/wsxu/Desktop/Projects/Stress/Toolkit/SetSteady'

#====get dir list
TList = sorted(glob.glob('../T*'), reverse=False)
TList = [Dir.split('/')[-1] for Dir in TList]

Start = TList.index('{}'.format(sys.argv[1]))
Final = TList.index('{}'.format(sys.argv[2]))

TList = [Dir for Dir in TList[Start:Final+1]]

#====tau data
for DirT in TList:
	call('cd ../{}/Prop/ && python {}/GetProp.py'.format(DirT, CODE_PATH), shell=True)

call('python GetData.py'.format(CODE_PATH), shell=True)
call('cd ../Figure/Dynamic/ && python Plot.py', shell=True)
call('cd ../Figure/Glass/ && python Plot.py', shell=True)
