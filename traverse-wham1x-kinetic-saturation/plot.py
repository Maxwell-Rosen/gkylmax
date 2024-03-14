import numpy as np
import matplotlib.pyplot as plt

# Read the data from the text file into a numpy array
Nz_vals_8 = np.loadtxt('Nz_vals_8.txt')
execution_times_8 = np.loadtxt('execution_times_8.txt')

Nz_vals_20 = np.loadtxt('Nz_vals_20.txt')
execution_times_20 = np.loadtxt('execution_times_20.txt')

Nz_vals_48 = np.loadtxt('Nz_vals_48.txt')
execution_times_48 = np.loadtxt('execution_times_48.txt')

# Run scanning Nz
# 5.24328 9.33268 13.0298 17.7237 24.3106 30.8948 38.3161 37.5154 47.9875 56.6498 64.6352 72.0228 82.9414 90.3513
# 40 60 80 100 120 140 160 180 200 220 240 260 280 300

# Run scanning vpar
print("8 GPUs")
print(Nz_vals_8)
print(execution_times_8)

print("20 GPUs")
print(Nz_vals_20)
print(execution_times_20)

print("48 GPUs")
print(Nz_vals_48)
print(execution_times_48)

plt.plot(Nz_vals_8 / 8.0, execution_times_8, 'bo-')
plt.plot(Nz_vals_20 / 20.0, execution_times_20, 'ro-')
plt.plot(Nz_vals_48 / 48.0, execution_times_48, 'go-')
plt.xlabel('Nz per GPU')
plt.ylabel('Execution time (s)')
plt.title('GPU saturation of WHAM1x')
plt.legend(['8 GPUs', '20 GPUs', '48 GPUs'])
plt.show()
