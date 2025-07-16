import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('phases_data.txt', delimiter=',', dtype=float)

cycle = data[:, 0]
print("Cycle data:", cycle)
ii = data[:, 1]
ie = data[:, 2]
tstart = data[:, 3]
tfinal = data[:, 4]
prev_frames = data[:, 5]
cycle_t_elc = data[:, 6]
cycle_t_ion = data[:, 7]

plt.figure(figsize=(10, 6))
for i in range(len(cycle)):
    colors = ['b', 'g', 'r', 'c']  # Assign colors for cycles 0,1,2,3
    if (cycle[i] == 0):
      plt.plot([tstart[i], tfinal[i]], [0, 0], marker='.', linestyle='-', color=colors[0])
    if (cycle[i] == 1):
      plt.plot([tstart[i], tfinal[i]], [1, 1], marker='.', linestyle='-', color=colors[1])
      plt.plot([tstart[i], tfinal[i]], [2, 2], marker='.', linestyle='-', color=colors[2])
    if (cycle[i] == 2):
      plt.plot([tstart[i], tfinal[i]], [1, 1], marker='.', linestyle='-', color=colors[1])
      plt.plot([tstart[i], tfinal[i]], [3, 3], marker='.', linestyle='-', color=colors[3])
    if (cycle[i] == 3):
      plt.plot([tstart[i], tfinal[i]], [3, 3], marker='.', linestyle='-', color=colors[3])

# Put text on the y axis for 0 as Electron RDP, 1 as Electron OAP, 2 as Ion RDP, and 3 as Ion OAP
plt.yticks([0, 1, 2, 3], ['Electron RDP', 'Electron OAP', 'Ion RDP', 'Ion OAP'])
# plt.grid(True)

plt.title('Phases Plot')
plt.xlabel('Time, μs')
plt.ylabel('Cycle type')
plt.show()