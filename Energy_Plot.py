import numpy as np
import matplotlib.pyplot as plt

data_file = np.loadtxt('test_v2_E_list.txt',delimiter=None)

time = data_file[:,0]

kine = data_file[:,1]

pote = data_file[:,2]

totE = data_file[:,3]

plt.plot(time,kine, label = "Kinetic Energy")
plt.plot(time, pote, label = "Potential Energy")
plt.plot(time, totE, label = "Total Energy", linestyle = '--', color ='r')
plt.legend()
plt.show()
plt.savefig('Energy_Plots.png')