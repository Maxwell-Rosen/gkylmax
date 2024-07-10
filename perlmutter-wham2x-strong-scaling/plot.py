import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def model(x, s):
    return 1 / (s + (1 - s) / x)

# Read the data from the text file into a numpy array
# num_procs = np.loadtxt('num_procs.txt')
# execution_times = np.loadtxt('execution_times.txt')

execution_times = np.array([ .1638, .0835, .0440, 0.0228, 0.0152, .0076, .0057])
num_procs = np.array([1, 2, 4, 8, 13, 52, 104])

speedup = execution_times[0]/execution_times

initial_guess = 0.05  # Initial guess for the parameter s
params, covariance = curve_fit(model, num_procs, speedup, p0=initial_guess)
s_optimal = params[0]
print(f"Optimal parameter s: {s_optimal}")


plt.plot(num_procs, speedup, 'bo', label='WHAM2x')
x_fit = np.linspace(min(num_procs), max(num_procs), 100)
y_fit = model(x_fit, s_optimal)
plt.plot(x_fit, y_fit, 'g-', label='Fit to WHAM2x')
y_ideal = model(num_procs, 0)
plt.plot(num_procs, y_ideal, 'r-', label='Ideal scaling')

plt.xlabel('Number of processes')
plt.ylabel('Fractional speedup relative to 1 proces')
# plt.title('Strong scaling of WHAM2x. GPU fraction is (1-s)=', s_optimal)
plt.title('Strong scaling of WHAM2x\n The code is ' + f"{1-s_optimal:.4f}" + ' fraction paralellized')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.savefig('wham2x-strong-scaling.png')
plt.close()

# Wall clock per time step
plt.plot(num_procs, execution_times, 'bo', label='WHAM2x')
plt.plot(x_fit, execution_times[0]/x_fit, 'r-', label='Ideal scaling')
plt.xlabel('Number of processes')
plt.ylabel('Wall clock time per time step (s)')
plt.xscale('log')
plt.yscale('log')
plt.title('Wall clock time per time step of WHAM2x')
plt.legend()
plt.savefig('wham2x-wall-clock.png')

