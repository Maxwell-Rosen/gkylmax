import numpy as np
import matplotlib.pyplot as plt

# Read the data from the text file into a numpy array
Nz_vals = np.loadtxt('Nz_vals.txt')
execution_times = np.loadtxt('execution_times.txt')

# Run scanning Nz
# 5.24328 9.33268 13.0298 17.7237 24.3106 30.8948 38.3161 37.5154 47.9875 56.6498 64.6352 72.0228 82.9414 90.3513
# 40 60 80 100 120 140 160 180 200 220 240 260 280 300

# Run scanning vpar
print(Nz_vals)
print(execution_times)

plt.plot(Nz_vals, execution_times, 'o-')
plt.xlabel('Nz')
plt.ylabel('Execution time (s)')
plt.title('Weak scaling of WHAM1x')
plt.show()
