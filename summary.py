import numpy as np
import matplotlib.pyplot as plt

file = open('results/mon_ver07/summary_problems.txt','w')

numbers = np.zeros(600)
offsets = np.zeros((600,30))
resolutions = np.zeros((600,30))
runlist = np.loadtxt("runs.dat",dtype=int)
i = -1

for run in runlist:
    if np.loadtxt("results/mon_ver07/st_time_res_"+str(run)+".txt").size>0:
        i += 1
        numbers[i] = run
        offsets[i,:] = np.loadtxt("results/mon_ver07/st_time_res_"+str(run)+".txt")[:,2]*1000
        resolutions[i,:] = np.loadtxt("results/mon_ver07/st_time_res_"+str(run)+".txt")[:,3]*1000
        if np.average(abs(offsets[i,:])) > 100:
            file.write(str(run)+" average offset is "+str("%.2f" % np.average(abs(offsets[i,:])))+" ps, greater than 100 ps\n")
            continue
        if np.average(resolutions[i,:]) > 400:
            file.write(str(run)+" average resolution is "+str("%.2f" % np.average(abs(resolutions[i,:])))+" ps, greater than 350 ps\n")
            continue
        if np.average(resolutions[i,:]) < 200:
            file.write(str(run)+" average resolution is "+str("%.2f" % np.average(abs(resolutions[i,:])))+" ps, smaller than 250 ps\n")
            continue
    else:
        file.write(str(run)+" file doesn't exist\n")

plt.scatter(numbers[:i],np.average(offsets[:i],axis=1))
plt.ylim([-200,200])
plt.xlabel("Run")
plt.ylabel("offset(ps)")
plt.title("offsets, average of all sectors")
plt.savefig('results/mon_ver07/summary_offsets_average.png')
plt.close()

fig, axs = plt.subplots(3,10,sharex=True,sharey=True,figsize=(30,10))
for j in range(3):
    axs[j,0].set_ylabel("offset(ps)")
    for k in range(10):
	axs[2,k].set_xlabel("Run")
	axs[j,k].scatter(numbers[:i],offsets[:i,10*j+k])
plt.ylim([-200,200])
plt.suptitle("offsets, individual sectors",fontsize=20)
plt.savefig('results/mon_ver07/summary_offsets_sector.png')
plt.close()

plt.scatter(numbers[:i],np.average(resolutions[:i],axis=1))
plt.ylim([200,500])
plt.xlabel("Run")
plt.ylabel("resolution(ps)")
plt.title("resolutions, average of all sectors")
plt.savefig('results/mon_ver07/summary_resolutions_average.png')
plt.close()

fig, axs = plt.subplots(3,10,sharex=True,sharey=True,figsize=(30,10))
for j in range(3):
    axs[j,0].set_ylabel("resolution(ps)")
    for k in range(10):
	axs[2,k].set_xlabel("Run")
	axs[j,k].scatter(numbers[:i],resolutions[:i,10*j+k])
plt.ylim([200,500])
plt.suptitle("resolutions, individual sectors",fontsize=20)
plt.savefig('results/mon_ver07/summary_resolutions_sector.png')
plt.close()
