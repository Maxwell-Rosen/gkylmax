import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img 

m_i_div_m_e = 1836

K_wall = 10
K = 10**(np.linspace(0, np.log10(K_wall), 100))
ephi_div_Te_K10 = -(m_i_div_m_e**(1/3)*(1/K**(2/3) - 1/K_wall**(2/3)))
ephi_div_Te_K10 = ephi_div_Te_K10 - ephi_div_Te_K10[0]

K_wall300 = 300
K300 = 10**(np.linspace(0,np.log10(K_wall300), 100))
ephi_div_Te_K300 = -(m_i_div_m_e**(1/3)*(1/K300**(2/3) - 1/K_wall300**(2/3)))
ephi_div_Te_K300 = (ephi_div_Te_K300 - ephi_div_Te_K300[0])*1

test_expression = 1/2*np.log(K300)
test_expression = (test_expression - test_expression[0])

phi_sheath = 1/(K_wall**(2/3))*(m_i_div_m_e)**(1/3)

ryutov_10 = np.loadtxt('ephi_mirDTe-ephiDte.csv', delimiter=',')
ryutov_300 = np.loadtxt('Ryutov_2005-fig2-K300-ephi_mirDTe-ephiDte.csv', delimiter=',') 

plt.figure()
plt.plot(K, ephi_div_Te_K10, 'b', label = "Ryutov equations K=10")
plt.plot(K300, ephi_div_Te_K300, 'm--', label = "Ryutov equations K=300")
# plt.plot(K300, test_expression, 'k--', label = "Test expression")
plt.plot([K_wall, K_wall], [ephi_div_Te_K10[-1], ephi_div_Te_K10[-1] + phi_sheath], 'r--', label = "Sheath jump")
plt.plot(ryutov_10[:, 0], ryutov_10[:, 1], 'r', label = "Ryutov-2005 K=10")
plt.plot(ryutov_300[:, 0], ryutov_300[:, 1], 'g', label = "Ryutov-2005 K=300")
plt.plot(K300, 1/2*np.log(K300), 'k--', label = "ln(K)/2")
plt.xlabel(r'$B_{\rm {mir}}/B$')
plt.ylabel(r'$e(\phi_{mir} - \phi)/T_e$')
plt.xscale('log')
plt.yscale('linear')
plt.xlim([1, 1000])
# plt.ylim([0, 4.4])
plt.legend()
plt.savefig('figure.png')
plt.show()

# # plot the image Ryutov-2005-figure-2.png
# plt.figure()
# image = img.imread('Ryutov-2005-figure-2.png')
# plt.imshow(image, extent=[1, 1000, 0, 3.5], aspect='auto')
# # Hide the x and y ticks
# plt.xticks([])
# plt.yticks([])
# plt.savefig('temp_fig2.png')
# plt.close()

# # Plot temp_fig1 on top of temp_fig2
# plt.figure()
# image = img.imread('temp_fig2.png')
# xshift = -96
# yshift = -0.348
# plt.imshow(image, extent=[1+xshift, 1000+xshift, 0+yshift, 3.5+yshift], aspect='auto')
# image = img.imread('temp_fig1.png')
# plt.imshow(image, extent=[1, 1000, 0, 3.5], aspect='auto', alpha=0.5)
# # plt.show()