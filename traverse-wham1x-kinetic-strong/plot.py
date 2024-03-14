import numpy as np
import matplotlib.pyplot as plt

# Read the data from the text file into a numpy array
num_procs = np.loadtxt('num_procs_240.txt')
execution_times = np.loadtxt('execution_times_240.txt')


# Run scanning Nz
# 5.24328 9.33268 13.0298 17.7237 24.3106 30.8948 38.3161 37.5154 47.9875 56.6498 64.6352 72.0228 82.9414 90.3513
# 40 60 80 100 120 140 160 180 200 220 240 260 280 300

plt.plot(num_procs, execution_times, 'bo-', label='WHAM1x')
plt.plot(num_procs, execution_times[0]/num_procs, 'k--', label='ideal')
# set x and y  axies to log scale
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of processes')
plt.ylabel('Execution time per update (s)')
plt.title('Strong scaling of WHAM1x for 240 cells in z')
plt.legend()
plt.show()
