import numpy as np
import matplotlib.pyplot as plt
import sys

def disp_usage():
	print("Usage: rf_kappa.py [filename] [-options]")
	print("Options:")
	print("-ifw: Specify input filename containing the forward time trajectory")
	print("-ibw: Specify input filename containing the backward time trajectory")
	print("-iv: Specify input filename containing velocities")
	
def pbc_dx(Dx,box_len):
    dx=np.zeros(np.shape(Dx))
    traj_ids=np.arange(np.shape(Dx)[1])
    for traj_id in traj_ids:
        dx[:,traj_id] = np.array([min(delx,box_len-delx) for delx in Dx[:,traj_id]])

    return(dx)
    
def cal_dx(x,box_len):
    n_mols=np.shape(x)[1]
    #print(n_mols)
    Dx=np.array([abs(x[:,i+1]-x[:,i]) for i in range(0,n_mols,2)],dtype='float16')
    Dx=Dx.T
    dx=pbc_dx(Dx,box_len)
    return(dx)

def cal_vrc(v):
    n_mols=np.shape(v)[1]
    vq=np.array([v[:,i+1]-v[:,i] for i in range(0,n_mols,2)],dtype='float16')
    vq=vq.T
    return(vq)
 
def indicator_upd(traj,frame_id,axis,xr,xp):
	delxr=abs(traj[frame_id,axis,:]-xr)
	delxp=abs(traj[frame_id,axis,:]-xp)
	return(delxr,delxp)

def absorb_at_boundary(traj,xr,xp,width):
    import matplotlib.pyplot as poltu
    frame_ids=np.arange(0,np.shape(traj)[0])
    traj_counts=np.zeros((n_frames,4))
    kappa_t=np.zeros(n_frames)
    kappa_2t=np.zeros(n_frames)
    #print(frame_ids)

    for frame_id in frame_ids:
      print("PROCESSING FRAME NUMBER:",frame_id,end="\r")
      delxrf,delxpf=indicator_upd(traj,frame_id,0,xr,xp)
      delxrb,delxpb=indicator_upd(traj,frame_id,1,xr,xp)

      prod_f=np.where(delxpf<=width)
      rec_f=np.where(delxrf<=width) 
      prod_b=np.where(delxpb<=width)
      rec_b=np.where(delxrb<=width)

      traj[frame_id+1:,0,prod_f[0]]=traj[frame_id,0,prod_f[0]]
      traj[frame_id+1:,1,prod_b[0]]=traj[frame_id,1,prod_b[0]]
      traj[frame_id+1:,0,rec_f[0]]=traj[frame_id,0,rec_f[0]]
      traj[frame_id+1:,1,rec_b[0]]=traj[frame_id,1,rec_b[0]]
    #print(traj)
    print("COMPLETED")	
    #print(prod_f,rec_f,prod_b,rec_b)
    #print(np.shape(prod_f))
    #plt.plot(frames,traj[:,0,prod_f[0][3]])  
    #plt.plot(frames,traj[:,0,:])
    #plt.show()
    
def track_traj(traj,frame_id,axis,xr,xp):
    delxr=traj[frame_id,axis,:]-xr
    delxp=traj[frame_id,axis,:]-xp
    return(delxr,delxp)
    
def cal_tc_gen(vel,rp,pr,rr,pp,print_flag):
    R=8.3143
    temp=300
    prob=np.zeros(np.shape(vel)[0])
    #vel=velocities[0,:]
    #prob=np.exp(-potential(traj[0,0,:],V0,n,L)/(R*temp))
    #prob=np.exp(-pot/(R*temp))
    prob[:]=1
    q=np.zeros(np.shape(vel)[0])
    q[rp]=1;q[pr]=-1;q[rr]=0;q[pp]=0
    #kplus=np.sum(vel[rp]*prob[rp]*q[rp])/(np.sum(abs(vel[rp])*prob[rp])+np.sum(abs(vel[rr])*prob[rr]))
    #kminus=np.sum(vel[pr]*prob[pr]*q[pr])/(np.sum(abs(vel[pr])*prob[pr])+np.sum(abs(vel[pp])*prob[pp]))
    flux_num=np.sum(vel*prob*q)
    flux_denom=np.sum(abs(vel)*prob)
    kappa=flux_num/flux_denom
    if(print_flag):
        print(np.shape(rp)[0],np.shape(pr)[0],np.shape(rr)[0],np.shape(pp)[0])
    #print(flux_num,flux_denom,kappa)
    #print(kplus,kminus,kplus+kminus,kplus/kminus,kminus/)
    return(kappa)#,kplus,kminus)

def analyzer_gen(n_frames,traj,vel,xr,xts,xp,width,stride,plot_flag,print_flag,hist_flag):
    traj_counts=np.zeros((n_frames,4))
    kappa=np.zeros(n_frames)
    kappa_plus=np.zeros(n_frames)
    kappa_minus=np.zeros(n_frames)
    traj_tot=np.zeros(n_frames)
    frame_ids=np.arange(0,np.shape(traj)[0])

    for frame_id in frame_ids:
      delxrf,delxpf=track_traj(traj,frame_id,0,xr,xp)
      delxrb,delxpb=track_traj(traj,frame_id,1,xr,xp)
      
      if(xp>xr):  
          xf=traj[frame_id,0,:]
          xb=traj[frame_id,1,:]
      else:
          xf=traj[frame_id,1,:]
          xb=traj[frame_id,0,:]

        
      rp=np.intersect1d(np.where(delxpf<=width),np.where(delxrb<=width)) #COUNTING R-->P TRAJECTORIES
      pr=np.intersect1d(np.where(delxpb<=width),np.where(delxrf<=width))#COUNTING P-->R TRAJECTORIES
      rr=np.intersect1d(np.where(delxrf<=width),np.where(delxrb<=width))#COUNTING R-->R TRAJECTORIES
      pp=np.intersect1d(np.where(delxpf<=width),np.where(delxpb<=width))#COUNTING P-->P TRAJECTORIES
        
      #rp=np.intersect1d(np.where(xf>=xp-width),np.where(xb<=xr+width)) #COUNTING R-->P TRAJECTORIES
      #pr=np.intersect1d(np.where(xb>=xp-width),np.where(xf<=xr+width))#COUNTING P-->R TRAJECTORIES
      #rr=np.intersect1d(np.where(xb<=xr+width),np.where(xf<=xr+width))#COUNTING R-->R TRAJECTORIES
      #pp=np.intersect1d(np.where(xb>=xp-width),np.where(xf>=xp+width))#COUNTING P-->P TRAJECTORIES

      #rp=np.intersect1d(np.where(xf>traj[0,0,:]),np.where(xb<traj[0,1,:])) #COUNTING R-->P TRAJECTORIES
      #pr=np.intersect1d(np.where(xb>traj[0,1,:]),np.where(xf<traj[0,0,:]))#COUNTING P-->R TRAJECTORIES
      #rr=np.intersect1d(np.where(xb<traj[0,1,:]),np.where(xf<traj[0,0,:]))#COUNTING R-->R TRAJECTORIES
      #pp=np.intersect1d(np.where(xb>traj[0,1,:]),np.where(xf>traj[0,0,:]))#COUNTING P-->P TRAJECTORIES
      
      rp=np.intersect1d(np.where(xf>xts),np.where(xb<xts)) #COUNTING R-->P TRAJECTORIES
      pr=np.intersect1d(np.where(xb>xts),np.where(xf<xts))#COUNTING P-->R TRAJECTORIES
      rr=np.intersect1d(np.where(xb<xts),np.where(xf<xts))#COUNTING R-->R TRAJECTORIES
      pp=np.intersect1d(np.where(xb>xts),np.where(xf>xts))#COUNTING P-->P TRAJECTORIES
    
      traj_counts[frame_id,0]=np.shape(rp)[0]
      traj_counts[frame_id,1]=np.shape(pr)[0]
      traj_counts[frame_id,2]=np.shape(rr)[0]
      traj_counts[frame_id,3]=np.shape(pp)[0]


      traj_tot[frame_id]=np.sum(traj_counts[frame_id,:])  
      #show_traj(rp)
      #traj_tot=traj_counts[frame_id,0]+traj_counts[frame_id,1]+traj_counts[frame_id,2]+traj_counts[frame_id,3]
      #if(traj_tot[frame_id]>0):
      kappa[frame_id]=cal_tc_gen(vel,rp,pr,rr,pp,print_flag)
      #kappa[frame_id],kappa_plus[frame_id],kappa_minus[frame_id]=cal_tc_gen(vel,rp,pr,rr,pp,print_flag)
      
      
      #print(traj_tot[frame_id],traj_counts[frame_id,0],traj_counts[frame_id,1],traj_counts[frame_id,2],traj_counts[frame_id,3],kappa[frame_id])
    if(plot_flag==True):
        #plt.xlim(180,300)  
        #plt.ylim(-10,25)
        plt.xlabel('Time(ps)')
        plt.ylabel('Trajectory Counts')
        plt.plot(frame_ids/stride,traj_counts[:,0],label='RP')
        plt.plot(frame_ids/stride,traj_counts[:,1],label='PR')
        plt.plot(frame_ids/stride,traj_counts[:,2],label='RR')
        plt.plot(frame_ids/stride,traj_counts[:,3],label='PP')
        #plt.plot(frame_ids/stride,traj_counts[:,3]+traj_counts[:,2],label='RR+PP')
        #plt.plot(frame_ids,traj_tot,label='Total')
        #plt.plot(frame_ids,traj_counts[:,0]-traj_counts[:,1],label='RP-PR')
        plt.legend()
        plt.grid()
        plt.show()  
        import matplotlib.pyplot as drw
        #drw.ylim(-0.0,0.0)
        #drw.xlim(180,300)
        drw.xlabel('Time(ps)')
        drw.ylabel('Transmission Coefficient ($\kappa$)')
        drw.plot(frame_ids/stride,kappa,label='$\kappa$')
        #drw.plot(frame_ids/stride,kappa_plus,label='$\kappa_{+}$')
        #drw.plot(frame_ids/stride,kappa_minus,label='$\kappa_{-}$')
        drw.grid()
        drw.legend()
        drw.savefig('kappa_vs_t.png',dpi=300)
        drw.show()
    if(hist_flag==True):   
        rx=np.concatenate((rr,pp))
        print(rx)
        trj,tp1=turning_points(xr,traj[:,:,rr])
        trj,tp2=turning_points(xp,traj[:,:,pp])
        tp=np.concatenate((tp1,tp2))
        print(tp)
        #plt.xlim(0,20)
        plt.xlabel('x ($A^{o}$)')
        plt.ylabel('No of Recrossings')
        plt.hist(tp,100)
        plt.savefig('Recrossing_Dist_Hist.png',dpi=300)
        plt.show()
    
    print(traj_counts[frame_id,0],traj_counts[frame_id,1],traj_counts[frame_id,2],traj_counts[frame_id,3],kappa[frame_id])
    return(kappa[frame_id],rp,pr,rr,pp)
    

n_args=len(sys.argv)
scale=0.5
stride=0.0001
ifname="force.fx"
ofname="friction_kernel.dat"
box_len=100
pos_rw=1.5
pos_pw=3.0

if(n_args>1):
	for i in range(1,n_args):
		if(sys.argv[i]=="-ifw"):
			ifwname=str(sys.argv[i+1])
		if(sys.argv[i]=="-iv"):
			ivname=str(sys.argv[i+1])
		if(sys.argv[i]=="-ibw"):
			ibwname=str(sys.argv[i+1])
		
else:
		disp_usage()
		sys.exit()


traj_fw=np.loadtxt(str(ifwname))
traj_bw=np.loadtxt(str(ibwname))

n_frames=np.shape(traj_fw)[0]
n_mols=np.shape(traj_fw)[1]
n_traj=n_mols

print(n_frames,n_traj)
frames=np.arange(n_frames)
traj_ids=np.arange(n_traj)

#traj_fw=cal_dx(traj_fw,box_len)
#traj_bw=cal_dx(traj_bw,box_len)

print(np.shape(traj_fw))
print(np.shape(traj_bw))

velocities=np.loadtxt(str(ivname),dtype='float16')
vq=velocities 

traj=np.zeros((n_frames,2,n_traj))
traj[:,0,:]=traj_fw[:,:]
traj[:,1,:]=traj_bw[:,:]

absorb_at_boundary(traj,pos_rw,pos_pw,0.05)
kappa,rp,pr,rr,pp=analyzer_gen(n_frames,traj,vq[0,:],pos_rw,2.25,pos_pw,0,10000,False,False,False)       
    
file=open('kappa.out','a')
file.write(str(kappa))
file.write("\n")
file.close()
