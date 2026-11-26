import os, sys
import numpy as np
from subprocess import call
from scipy.optimize import curve_fit
from scipy.signal import argrelmax
from scipy.interpolate import interp1d

Job = os.getcwd().split('/')[-1]

#====function Get_Single: process single-line properties
def Get_Single(DirList, NDir, Fname):
#	====get information about the data structure
	Data = np.loadtxt('../{}/{}/{}'.format(DirList[0], Job, Fname))
	N = Data.size

#	====combine data from different files
	Prop = np.zeros((N, NDir))
	for i in range(len(DirList)):
		Data = np.loadtxt('../{}/{}/{}'.format(DirList[i], Job, Fname))
		Prop[:,i] = Data

#	====get average and standard deviation	
	Mean, STD = np.mean(Prop, axis=1), np.std(Prop, axis=1)

#	====put the result in an array and save it	
	Result = np.zeros((1, 2*N))
	for n in range(N):
		Result[0,2*n+0] = Mean[n]
		Result[0,2*n+1] = STD[n]	
	np.savetxt(Fname, Result, fmt='%17.8e')
		
#	print('{:>10} - File & Item: {:<6}  {:>1}'.format('Single', Fname, N))
	
#====function Get_Static: process static properties
def Get_Static(DirList, NDir, Fname):
#	====get information about the data structure
	Data = np.loadtxt('../{}/{}/{}'.format(DirList[0], Job, Fname))
	X = Data[:,0]
	NLine = Data.shape[0]
	NProp = Data.shape[1] - 1
	
#	====combine data from different files
	Prop = np.zeros((NLine, NProp, NDir))
	for i in range(len(DirList)):
		Data = np.loadtxt('../{}/{}/{}'.format(DirList[i], Job, Fname))
		for n in range(NProp):
			Prop[:,n,i] = Data[:,n+1]

#	====get average and standard deviation
	Mean = np.zeros((NLine, NProp))
	STD = np.zeros((NLine, NProp))
	for n in range(NProp):
		Mean[:,n], STD[:,n] = np.mean(Prop[:,n,:], axis=1), np.std(Prop[:,n,:], axis=1)

#	====put the result in an array and save it		
	Result = np.zeros((NLine, 2*NProp+1))
	Result[:,0] = X
	for n in range(NProp):
		Result[:,2*n+1] = Mean[:,n]
		Result[:,2*n+2] = STD[:,n]	
	np.savetxt(Fname, Result, fmt='%17.8e')
	
#	print('{:>10} - File & Item: {:<6}  {:>1}'.format('Static', Fname, NProp))
	
#====function Get_Dynamic: process dynamic properties
def Get_Dynamic(DirList, NDir, Fname, Del):
#	====put data in one file
	for i in range(len(DirList)):
		call('>{}'.format(DirList[i]), shell=True)
		Xmax = 1e9
		for Sub in ['L', 'M', 'S']:
			Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[i], Job, Fname, Sub))
			Part = Data[Data[:,0]<Xmax]
			Min = 1.0
			Max = Part.shape[0]
			Num = np.log10(Max/Min)/Del
			Num = Num.astype(int)
			Index = np.unique(np.geomspace(Min, Max, num=Num, dtype=int)) - 1
			Part = Part[Index]	
			np.savetxt('tmp', Part, fmt='%17.8e')
			call('cat tmp >>{}'.format(DirList[i]), shell=True)
			Xmax = 1.0*(Data[1,0] - Data[0,0])	
#		====remove duplicates and sort data
		Data = np.loadtxt(DirList[i])
		Sort = Data[np.argsort(Data[:,0], kind='heapsort')]
		Value, Index = np.unique(Sort[:,0], return_index=True)
		Data = Sort[Index]
		np.savetxt(DirList[i], Data, fmt='%17.8e')
	
#	====get information about the data structure
	Data = np.loadtxt(DirList[0])
	X = Data[:,0]
	NLine = Data.shape[0]
	NProp = Data.shape[1] - 1
	
#	====combine data from different files
	Prop = np.zeros((NLine, NProp, NDir))
	for i in range(len(DirList)):
		Data = np.loadtxt(DirList[i])
		for n in range(NProp):
			Prop[:,n,i] = Data[:,n+1]

#	====get average and standard deviation
	Mean = np.zeros((NLine, NProp))
	STD = np.zeros((NLine, NProp))
	for n in range(NProp):
		Mean[:,n], STD[:,n] = np.mean(Prop[:,n,:], axis=1), np.std(Prop[:,n,:], axis=1)

#	====put the result in an array and save it		
	Result = np.zeros((NLine, 2*NProp+1))
	Result[:,0] = X
	for n in range(NProp):
		Result[:,2*n+1] = Mean[:,n]
		Result[:,2*n+2] = STD[:,n]	
	np.savetxt(Fname, Result, fmt='%17.8e')
	
#	====remove tmp files
	call('rm tmp', shell=True)
	for i in range(len(DirList)):
		call('rm {}'.format(DirList[i]), shell=True)
	
#	print('{:>10} - File & Item: {:<6}  {:>1}'.format('Dynamic', Fname, NProp))

#====function Get_Tau: process Tau from FSQ
def Get_Tau(DirList, NDir, Fname, ID, FTau):
#	====number of items	
	N = 2 # 1st--position, 2nd--value
	
#	====combine data from different files
	Prop = np.zeros((N, NDir))
	for i in range(len(DirList)):
		k = 0
		FType = 'L'
		Prop[0,i] = 1e6
		Prop[1,i] = 1e0
		for Sub in ['S', 'M', 'L']:
			Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[i], Job, Fname, Sub))
			NLine = len(Data[:,0])
			Lower = NLine - 50
			k = np.argmin(np.abs(Data[:,ID]-FTau))
			Test = Data[Lower:NLine,ID]
			if k<Lower and (Test<FTau).all():
				X = Data[:,ID]
				Y = Data[:,0]
				Value, Index = np.unique(X, return_index=True)
				finter = interp1d(X[Index], Y[Index], kind='linear')
				Prop[0,i] = finter(FTau)
				Prop[1,i] = FTau
				FType = Sub
				break
		
#	====get average and standard deviation	
	Mean, STD = np.mean(Prop, axis=1), np.std(Prop, axis=1)

#	====put the result in an array and save it	
	Result = np.zeros((1, 2*N))
	for n in range(N):
		Result[0,2*n+0] = Mean[n]
		Result[0,2*n+1] = STD[n]	
	np.savetxt('Tau{}_{:03d}_{}'.format(ID, int(round(1e2*FTau)), Fname), Result, fmt='%17.8e')
	
#	print('{:>10} - File: {}_{}'.format('Tau{:03d}'.format(int(round(1e2*FTau))), Fname, FType))
	
	return (FType)

#====function Get_Peak: process peak from a dynamic property
def Get_Peak(DirList, NDir, Fname, ID, FType):	
#	====number of items	
	N = 2 # 1st--peak position, 2nd--peak magnitude
	
#	====combine data from different files
	Prop = np.zeros((N, NDir))
	for i in range(len(DirList)):
		Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[i], Job, Fname, FType))
		if Fname=='MSD':
			Index = sorted(np.unique(Data[:,ID], return_index=True)[1])
			List = argrelmax(Data[Index,ID], order=20)[0]
			if len(List)>0:
				if len(List)==1:
					k = List[0]
				elif len(List)==2:
					k = List[np.argmax(Data[List[0:2],ID])]
				else:
					k = List[np.argmax(Data[List[0:3],ID])]
			else:
				k = np.argmax(Data[:,ID])
		else:
			k = np.argmax(Data[:,ID])
		Prop[0,i] = Data[k,0]
		Prop[1,i] = Data[k,ID]
	
#	====get average and standard deviation	
	Mean, STD = np.mean(Prop, axis=1), np.std(Prop, axis=1)

#	====put the result in an array and save it	
	Result = np.zeros((1, 2*N))
	for n in range(N):
		Result[0,2*n+0] = Mean[n]
		Result[0,2*n+1] = STD[n]	
	np.savetxt('Peak{}_{}'.format(ID, Fname), Result, fmt='%17.8e')
	
#	print('{:>10} - File: {}_{}'.format('Peak', Fname, FType))
	
#====function Get_DWF: process Debye-Waller factor from MSD
def Get_DWF(DirList, NDir, Fname, ID):
#	====number of items	
	N = 2 # 1st--position, 2nd--value

#	====determine the file to be used
	for Sub in ['L', 'M', 'S']:
		Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[0], Job, Fname, Sub))
		X = Data[:,0]
		DelX = X[1] - X[0]
		if DelX<0.2:
			break
		
#	====get index corresponding to DWF
	k = np.argmin(np.abs(X - 2.0))
	
#	====combine data from different files
	Prop = np.zeros((2,NDir))
	for i in range(len(DirList)):
		Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[i], Job, Fname, Sub))
		Prop[0,i] = Data[k,0]
		Prop[1,i] = Data[k,ID]
		
#	====get average and standard deviation	
	Mean, STD = np.mean(Prop, axis=1), np.std(Prop, axis=1)

#	====put the result in an array and save it	
	Result = np.zeros((1, 2*N))
	for n in range(N):
		Result[0,2*n+0] = Mean[n]
		Result[0,2*n+1] = STD[n]	
	np.savetxt('DWF{}_{}'.format(ID,Fname), Result, fmt='%17.8e')
	
#	print('{:>10} - File: {}_{}'.format('DWF', Fname, Sub))
	
#====function Get_DR2: process derivative from MSD
def Get_DR2(DirList, NDir, Fname, ID, Del):
	Fout = 'DR2{}_{}'.format(ID,Fname)
	
#	====loop over different files
	call('>{}'.format(Fout), shell=True)
	Xmax = 1.0e9
	for Sub in ['L', 'M', 'S']:
		Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[0], Job, Fname, Sub), skiprows=1)
		X = Data[:,0]
		DelX = X[1] - X[0]
		NLine = Data.shape[0]
	
#		====combine data from different files
		Prop = np.zeros((NLine-2, NDir))
		for i in range(len(DirList)):
			Data = np.loadtxt('../{}/{}/{}_{}'.format(DirList[i], Job, Fname, Sub), skiprows=1)
			for j in range(1, NLine-1):
				Prop[j-1,i] = (Data[j,0] + 0.5*DelX)*np.log(Data[j+1,1]/Data[j,1])/DelX
		
#		====get average and standard deviation	
		Mean, STD = np.mean(Prop, axis=1), np.std(Prop, axis=1)

#		====put the result in an array
		Result = np.zeros((NLine-2, 3))
		Result[:,0] = X[1:NLine-1] + 0.5*DelX
		Result[:,1] = Mean[:]
		Result[:,2] = STD[:]
		Min = 1.0
		Max = Result.shape[0]
		Num = np.log10(Max/Min)/Del
		Num = Num.astype(int)
		Index = np.unique(np.geomspace(Min, Max, num=Num, dtype=int)) - 1
		Result = Result[Index]

#		====remove strange values	
		f = open('{}_{}'.format(Fout, Sub), 'w')
		for i in range(1, Result.shape[0]-1):
			y1 = (Result[i,1] - Result[i-1,1])/Result[i-1,1]
			y2 = (Result[i,1] - Result[i+1,1])/Result[i+1,1]
			if (np.abs(y1)<0.05 and np.abs(y2)<0.05) or y1*y2<0.0:
				if Result[i,0]<5.0*Xmax:
					f.write('{:17.8e} {:17.8e} {:17.8e}\n'.format(Result[i,0], Result[i,1], Result[i,2]))
		f.close()
		
		call('cat {}_{} >>{}'.format(Fout, Sub, Fout), shell=True)
		call('rm {}_{}'.format(Fout, Sub), shell=True)
		Xmax = 5.0*DelX

#	====save final result	
	Data = np.loadtxt(Fout)
	Sort = Data[np.argsort(Data[:,0], kind='heapsort')]
	Value, Index = np.unique(Sort[:,0], return_index=True)
	Data = Sort[Index]
	np.savetxt(Fout, Data, fmt='%17.8e')
	
#	print('{:>10} - File: {}'.format(Fout, Fname))

#====function Get_Hyper: process hyperuniform properties from SSF	
def Get_Hyper(DirList, NDir, Fname):
#	====fitting function
	def func(x, a, b):
		return (a + b*x)

#	====number of items	
	N = 5 # 1st--Qp, 2nd--Sp, 3rd--S0, 4th--B0, 5th--Sr
	
#	====combine data from different files
	Prop = np.zeros((N, NDir))
	for i in range(len(DirList)):
		Data = np.loadtxt('../{}/{}/{}'.format(DirList[i], Job, Fname))
		X = Data[:,0]
#		====get peak
		Part = Data[X>2.0]
		k = np.argmax(Part[:,1])
		Prop[0,i] = Part[k,0]
		Prop[1,i] = Part[k,1]
#		====fit at low Q
		Part = Data[X<1.0]
		popt, pcov = curve_fit(func, Part[:,0], Part[:,1])
		Prop[2,i] = popt[0]
		Prop[3,i] = popt[1]
#		====ratio
		Prop[4,i] = Prop[2,i]/Prop[1,i]
		
#	====get average and standard deviation	
	Mean, STD = np.mean(Prop, axis=1), np.std(Prop, axis=1)

#	====put the result in an array and save it	
	Result = np.zeros((1, 2*N))
	for n in range(N):
		Result[0,2*n+0] = Mean[n]
		Result[0,2*n+1] = STD[n]	
	np.savetxt('Hyper', Result, fmt='%17.8e')
	
#	print('{:>10} - File: {}'.format('Hyper', Fname))
	
#====function Get_LpLb: process persistence length from BCF
def Get_LpLb(DirList, NDir, Fname):
#	====fitting function
	def func(x, a):
		return np.exp(-x/a)
	
#	====combine data from different files
	Prop = np.zeros(NDir)
	for i in range(len(DirList)):
		Data = np.loadtxt('../{}/{}/{}'.format(DirList[i], Job, Fname))
		popt, pcov = curve_fit(func, Data[:,0], Data[:,1])
		Prop[i] = popt[0]
		
#	====get average and standard deviation	
	Mean, STD = np.mean(Prop), np.std(Prop)

#	====put the result in an array and save it	
	Result = np.zeros((1, 2))
	Result[0,0] = Mean
	Result[0,1] = STD	
	np.savetxt('LpLb', Result, fmt='%17.8e')
	
#	print('{:>10} - File: {}'.format('LpLb', Fname))
