from numpy import *
from MDAnalysis import *
import os,sys
import matplotlib.pyplot as plt
import matplotlib
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

def Dist(a, b):
    r = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
    return r


name = str(sys.argv[1]) # Name of your SMD output directory
sel_const = str(sys.argv[2]) # First selection (to calculate center of mass and H-bonds) - constrained part
sel_pull = str(sys.argv[3]) # Second selection (to calculate center of mass and H-bonds) - pulled part

extract_forces = False # change it to False after force files are generated
calculate_HB = False # change it to False after HB files are generated


list_of_dirs = [subdir for subdir in os.listdir(name) if subdir[0:9] == 'SMD_theta']

## this part takes some time, since it is analysing your SMD trajectories looking for possible H-bonds between selections
u = Universe(name+'/SMD_constraints.pdb')
fix = u.select_atoms(sel_const).center_of_mass() #
bunch = [] #here the information about pulling vectors will be stored

for dir in list_of_dirs:
    print(dir," calculations started.")
    f = open(name+'/'+dir+'/mdrun.log','r')


    # SMD force vs time and SMD force vs distance
    if extract_forces:

        ft = open(name + '/' + dir + '/smd_force_time.dat', 'w')
        ft.write('# time [fs] vs force of pulling\n')
        fd = open(name + '/' + dir + '/smd_force_dist.dat', 'w')
        fd.write('# distance between COMs vs force of pulling\n')

        for line in f.readlines():
            if len(line)>5:
                if line[:5] == 'SMD  ':
                    t = line.split()[1] #timestep
                    r = array([float(line.split()[2]),float(line.split()[3]),float(line.split()[4])]) # COM position of pulled domain
                    f = array([float(line.split()[5]),float(line.split()[6]),float(line.split()[7])]) # Pulling force, starting from COM
                    ft.write(t+' '+str(((f[0]-r[0])**2 + (f[1]-r[1])**2 + (f[2]-r[2])**2)**0.5)+'\n')
                    fd.write(str(Dist(r,fix))+' '+str(((f[0]-r[0])**2 + (f[1]-r[1])**2 + (f[2]-r[2])**2)**0.5)+'\n')
        ft.close()
        fd.close()

    if calculate_HB:
        ht = open(name + '/' + dir + '/smd_hb_time.dat', 'w')
        ht.write('# time [fs] vs number of hydrogen bonds between domains\n')
        # Hydrogen bonds number vs time
        psf = [file for file in os.listdir(name) if file[-3:] == 'psf']
        print('Hydrogen bonds calculating - that may take some time...')
        u = Universe(name+'/'+psf[0] ,name+'/'+dir+'/md.dcd')
        hbonds = HBA(universe=u,hydrogens_sel="protein and name H*",donors_sel=sel_const,acceptors_sel=sel_pull)
        hbonds.run()
        out = hbonds.count_by_time()
        hbonds_rev = HBA(universe=u,hydrogens_sel="protein and name H*",donors_sel=sel_pull,acceptors_sel=sel_const)
        hbonds_rev.run()
        out2 = hbonds_rev.count_by_time()
        t=0
        for i in (out+out2):
            ht.write(str(t)+' '+str(i)+'\n')
            t+= 0.05 #[ns]
        ht.close()



    #Writing down used pulling vectors
    theta = double(dir.split('_')[2])
    phi = double(dir.split('_')[4])
    x = cos(deg2rad(phi))*sin(deg2rad(theta))
    y = sin(deg2rad(phi))*sin(deg2rad(theta))
    z = cos(deg2rad(theta))
    bunch.append([x,y,z])
savetxt(name+'/bunch_of_vectors.dat',bunch)




# Plots - this script is adjusted for 9 pulling directions. If you want more, please increase the number of rows (9) and number of columns (3) accordingly
fig, ax = plt.subplots(9,3,sharex='col',figsize=(5,10))

i=0
sorted=['SMD_theta_0_phi_0','SMD_theta_45_phi_0','SMD_theta_45_phi_90','SMD_theta_45_phi_180','SMD_theta_45_phi_270','SMD_theta_90_phi_0','SMD_theta_90_phi_90','SMD_theta_90_phi_180','SMD_theta_90_phi_270'

]
# looking for a maximum force value and HB numbers in all files to rescale figures
max_force = max( [max(loadtxt(name+'/'+dir+'/smd_force_time.dat')[:,1]) for dir in sorted[:]] )
max_HB = max( [max(loadtxt(name+'/'+dir+'/smd_hb_time.dat')[:,1]) for dir in sorted[:]] )
print(max_force, max_HB)
for dir in sorted[:]:
    print(dir)

    ft = loadtxt(name+'/'+dir+'/smd_force_time.dat')
    ht = loadtxt(name+'/'+dir+'/smd_hb_time.dat')
    fd = loadtxt(name+'/'+dir+'/smd_force_dist.dat')

    #Force vs time
    ax[i,0].plot(ft[:,0]/500000,ft[:,1],linewidth=2)
    ax[i,0].set_xlim(0,)
    ax[i, 0].set_ylim(0, max_force)
    #ax[i,0].set_xlabel('Time [ns]').set_fontsize(16)
    #ax[i,0].set_xticks(14)
    #ax[i,0].set_ylabel('Force [pN]').set_fontsize(16)
    #ax[i,0].yticks(size=14)


    #Force vs distance
    ax[i,1].plot(fd[:,0],fd[:,1],'.')
    ax[i,1].set_xlim(30,55)
    ax[i,1].set_ylim(0, max_force)
    #ax[i,1].set_xlabel(r'Distance [$\AA$]').set_fontsize(16)

    #ax[i,1].set_ylabel('Force [pN]').set_fontsize(16)



    #hydrogen bonds number vs time
    ax[i,2].plot(ht[:,0]*2,ht[:,1],linewidth=2)
    ax[i,2].set_xlim(0,)
    ax[i,2].set_ylim(0,max_HB)
    #ax[i,2].set_xlabel('Time [ns]').set_fontsize(16)
    #ax[i,2].set_ylabel('HB number').set_fontsize(16)
    i+= 1
fig.tight_layout()
plt.savefig(name+'/all.png',dpi=600)


plt.show()
