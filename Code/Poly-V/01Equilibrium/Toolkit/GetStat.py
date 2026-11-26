import os, sys, glob
from GetFunc import *

Param = os.getcwd().split('/')[-2]

#====get dir list
DirList = sorted(glob.glob('../0*'))
DirList = [Dir.split('/')[-1] for Dir in DirList]
NDir = len(DirList)

print('{}  Sample number & List: {}  {}'.format(Param, NDir, DirList))

#====single-line properties
List = ['Thermo', 'Size', 'LB2', 'ANG']
for Fname in List:
	Get_Single(DirList, NDir, Fname)
	
#====static properties
List = ['RDF', 'SSF', 'BCF', 'MID', 'PLB2', 'PANG']
for Fname in List:
	Get_Static(DirList, NDir, Fname)

#====hyperuniform properties from SSF
Fname = 'SSF'
Get_Hyper(DirList, NDir, Fname)

#====persistence length from BCF
Fname = 'BCF'
Get_LpLb(DirList, NDir, Fname)
