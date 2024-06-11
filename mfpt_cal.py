#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

def disp_usage():
	print("Usage: mfpt_beta [filename] [-options]")
	print("Options:")
	print("-i: Specify matrix input filename")
	print("-icv: Specify COLVAR input filename")
	print("-o: Specify output filename")
	print("-rw: Specify reactant basin location")
	print("-pw: Specify product basin location")
	print("-scale: Specify relative viscosity")
	print("-stride: Specify stride at which frames were recorded")
	print("-dt: Specify simulation timestep(in fs)")
	print("-w: Specify basin width")
	print("-l: Specify box length for pbc correction")
	
	

n_args=len(sys.argv)
box_len=100
pos_rw=1.5
pos_pw=3.0

ofname='mfpt_vs_eta.dat'
scale=1.0
stride=1.0
dt=1.0
width=0.0
cv_flag=False

if(n_args>=2):
	for i in range(1,n_args):
		if(sys.argv[i]=="-i"):
			ifname=str(sys.argv[i+1])
		if(sys.argv[i]=="-icv"):
			cvfname=str(sys.argv[i+1])
			cv_flag=True	
		if(sys.argv[i]=="-o"):
			ofname=str(sys.argv[i+1])	
		if(sys.argv[i]=="-scale"):
			scale=float(sys.argv[i+1])
		if(sys.argv[i]=="-stride"):
			stride=float(sys.argv[i+1])
		if(sys.argv[i]=="-w"):
			width=float(sys.argv[i+1])
		if(sys.argv[i]=="-rw"):
			pos_rw=float(sys.argv[i+1])
		if(sys.argv[i]=="-pw"):
			pos_pw=float(sys.argv[i+1])
		if(sys.argv[i]=="-dt"):
			dt=float(sys.argv[i+1])
		if(sys.argv[i]=="-l"):
			box_len=float(sys.argv[i+1])		
else:
	disp_usage()
	sys.exit(0)

if(cv_flag):
	x=np.loadtxt(str(cvfname),comments="#")
	x=x[:,1:]	
else:
	x=np.loadtxt(str(ifname))


n_mols=np.shape(x)[1]


traj_ids=np.arange(np.shape(x)[1])
frames=np.arange(np.shape(x)[0])

print("No of trajectories:",np.shape(x)[1])
print("No of frames:",np.shape(x)[0])
print("Scale:", scale)
print("Stride:",stride)


traj=x

#MFPT CALCULATION
fpts=[]
for traj_id in traj_ids:
	diff=np.abs(traj[:,traj_id]-pos_pw)
	frnr=np.where(diff<=width)
	#print(frnr)
	if(len(frnr[0])>0):    
	    fpt=np.min(frnr[0][0])
	    fpts.append(fpt)

#print(fpts)	
mfpt=np.mean(fpts)
err=np.std(fpts)/np.sqrt(len(fpts))

MFPT=(mfpt*stride*dt)/1000
serr=(err*stride*dt)/1000

print("MFPT(ps):",MFPT)
print("Averaged over "+str(len(fpts))+" transitions")
#Writing to file
file=open(str(ofname),'a')
file.write(str(scale)+"\t"+str(MFPT)+"\t"+str(serr)+"\t"+str(len(fpts)))
file.write("\n")
file.close()

