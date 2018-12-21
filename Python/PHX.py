import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import MPC_Functions as mpcfunc

# Simulation and MPC Parameters
h = .1
Np = 10
Nc = 4
ref = 42
w_f = 10
N_sim = 20

# Continuous time matrices
ACT = np.matrix([[-1.4472, 0, 1.2904, 0],
                 [0, -1.3144, 1.2324, -.028897],
                 [2.9048, 2.9902, -5.895, 0],
                 [0, 0, 0, -.1]])
                
BwCT = np.matrix([[0],
                  [45.107],
                  [0],
                  [0]])
                 
BuCT = np.matrix([[0],
                  [0],
                  [0],
                  [.1]])
 
CCT = np.matrix([[8330000, 0, 0, 0]])   

# Discrete time matrices             
ADT = np.identity(4) + ACT*h + np.linalg.matrix_power(ACT,2)*h**2/2 + np.linalg.matrix_power(ACT,3)*h**3/6 \
                     +np.linalg.matrix_power(ACT,4)*h**4/24

BwDT = BwCT*h + np.dot(ACT, BwCT)*h**2/2 + 1/6*np.dot(np.linalg.matrix_power(ACT,2),BwCT)*h**3\
              + 1/24*np.dot(np.linalg.matrix_power(ACT,3),BwCT)*h**4

BuDT = BuCT*h + np.dot(ACT, BuCT)*h**2/2 + 1/6*np.dot(np.linalg.matrix_power(ACT,2),BuCT)*h**3\
              + 1/24*np.dot(np.linalg.matrix_power(ACT,3),BuCT)*h**4

CDT = CCT
              
# Calling function that augments matrices and carries out receding horizon control
#Phi_Phi, Phi_F, Phi_R, Ae, Be, Ce, k, u1, y1 = mpcfunc.mpcgain(ADT, BuDT, BwDT, CDT, Nc, Np, ref, w_f, N_sim)

E = np.matrix([[2, -1], [-1, 1]])
F = np.matrix([[-1], [0]])
M = np.matrix([[-2, 0], [0, -1], [3, 2]])
gamma = np.matrix([[0], [0], [4]])

x = mpcfunc.QPHild(E, F, M, gamma)


# Plot response
'''
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(k[0,:], y1[0,:])
ax2.plot(k[0,:], u1[0,:])
ax1.set_ylabel('Output')
ax2.set_xlabel('Sampling Instant')
ax2.set_ylabel('Control')
'''
