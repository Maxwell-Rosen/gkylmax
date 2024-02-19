import numpy as np
import matplotlib.pyplot as plt

# Read the data from the text file into a numpy array
Nz_vals = np.loadtxt('Nz_vals.txt')
execution_times = np.loadtxt('execution_times.txt')

print(Nz_vals)
print(execution_times)

plt.plot(Nz_vals, execution_times, 'o-')
plt.xlabel('Nz')
plt.ylabel('Execution time (s)')
plt.title('Weak scaling of WHAM1x')
plt.show()
