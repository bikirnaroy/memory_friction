#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

def disp_usage():
	print("Usage: facf_cal.py [filename] [-options]")
	print("Options:")
	print("-i: Specify input filename")
	print("-o: Specify output filename")
	print("-T: Specify the temperature in Kelvin")
	print("-m: Provide mass (in kg per mol) to Mass correct the facf.")
	print("-scale: maximum lag for acf as a fraction of total no of frames.Default value is 0.5")
	print("-stride: Timestep at which the forces were recorded in ps ")
	print("-xmax: Maximum lag value till which to plot facf. Set -1 for full range")
	print("-plot: Plot the normalized friction kernel.")

def mass_corr(mu_eff,temp):
	R=8.3143 #J K-1 mole-1
	T=temp #K
	denom=(mu_eff*R*T)
	conv_fac=(4.18*4.18)*100 #applicable for LAMMPS "real" unit style, where forces have units kcal mol-1 Angstrom -1.Change it for other unit styles or MD engines.
	alpha=conv_fac/denom #0.0333684
	return(alpha)

def facf_cal_fast(data,scale,plot=False):
    
    n_traj=np.shape(data)[1]
    n_frames=np.shape(data)[0]
    
    
    lag_max=int(scale*n_frames)
    delta=np.arange(lag_max)
    
    frames=np.arange(n_frames)
    
    forces=data    
    facf=np.zeros(lag_max+1)
    print("Starting Calculation...")
    
    facf[0]=np.mean(forces*forces)

    for lag in range(1,lag_max+1):
        progress=lag*100/lag_max
        if(progress%10==0):
        	print("Calculation Completed:",progress,"%",end="\r") 
        pdt=forces[:-lag,:]*forces[lag:,:]
        facf[lag]=np.mean(pdt)
        #print(facf[lag])
    if(plot):
    	plt.plot(delta,facf)
    	plt.show() 
           

    return(facf)
    

n_args=len(sys.argv)
scale=0.5
stride=0.0001
ifname="force.fx"
ofname="friction_kernel.dat"
xmax=-1
plot=False
temp=300 #Temperature in Kelvin
mass=12*1.0e-3 #Mass of atoms constituting the reactive mode in kg/mol
mu_eff=6*1.0e-3 #kg mol-1 value assuming reduced  mass of 2 Carbon atoms bonded together 

if(n_args>1):
	for i in range(1,n_args):
		if(sys.argv[i]=="-i"):
			ifname=str(sys.argv[i+1])
		if(sys.argv[i]=="-scale"):
			scale=float(sys.argv[i+1])
		if(sys.argv[i]=="-stride"):
			stride=float(sys.argv[i+1])
		if(sys.argv[i]=="-xmax"):
			xmax=float(sys.argv[i+1])
		if(sys.argv[i]=="-plot"):
			plot=True
		if(sys.argv[i]=="-o"):
			ofname=str(sys.argv[i+1])
		if(sys.argv[i]=="-m"):
			mu_eff=float(sys.argv[i+1])	
		if(sys.argv[i]=="-T"):
			temp=float(sys.argv[i+1])	
	
    	
else:
		disp_usage()
		sys.exit()
    	
print("Proceeding with default parameters for unspecified flags.")
print(mu_eff)
fx=np.loadtxt(str(ifname))
n_frames=np.shape(fx)[0]
n_atoms=np.shape(fx)[1]
#force=np.sqrt(fx*fx+fy*fy+fz*fz)
force=fx

 

#n_frames=np.shape(force)[0]
n_mols=int(n_atoms/2)
f_rc=force.T
#f_rc=np.array([abs(force[:,i+1]-force[:,i]) for i in range(0,n_atoms,2)])
#f_rc=np.array([(force[:,i+1]-force[:,i]) for i in range(0,n_atoms,2)])
#f_rc=f_rc*(mu_eff/mass)
data=f_rc.T    
frames=np.arange(scale*n_frames)[:]*stride
#print(data)
facf=facf_cal_fast(data,scale) #Just the force autocorrelation function


print(mass_corr(mu_eff,temp))
facf*=mass_corr(mu_eff,temp) #friction kernel in ps-2 units

data2write=np.c_[frames,facf]
np.savetxt(ofname,data2write)

if(plot):
	if (xmax!=-1):
		plt.xlim(0,xmax)
	plt.xlabel('Time(ps)')
	plt.ylabel('$\zeta(t)(ps^{-2})$')
	plt.plot(frames,facf[:],'-')
	plt.savefig('facf.png',dpi=300)
	plt.show() 


