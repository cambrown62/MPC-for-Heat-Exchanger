
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
import MPC_Functions as mpcfunc
#import control
#from scipy.integrate import odeint

# System parameters and continuous time state matrices
zeta = .8
wn = 5.0
h = .01
Np = 10
Nc = 4
ref = 5
w_f = 20

ACT = np.matrix([[0., 1.],[-wn**2, -2*zeta*wn]])
BCT = np.matrix([[0.],[wn**2]])
CCT = np.matrix([[1., 0.]])
DCT = np.matrix([[0.]])

'''
# Continuous time input/output matrices for feedback linearization
AioCT = np.matrix([[0, 1],[0, 0]])
BioCT = np.matrix([[0],[1]])
CioCT = np.matrix([[1, 0]])

# Continuous time reference signal matrices
p1 = 1
p2 = 5

ArCT = np.matrix([[0, 1], [-p1*p2, -p1-p2]])
BrCT = np.matrix([[0],[p1*p2]])
'''

# Discrete time state matrices found using R-K in Matlab
ADT = np.matrix([[.99878, .0096065], [-.24016, .92193]])
BDT = np.matrix([[.0012171],[.24016]])
CDT = CCT
DDT = DCT

'''
# Discrete time input/output matrices
AioDT = np.matrix([[1, 0.01],[0, 1]])
BioDT = np.matrix([[5e-5],[.01]])
CioDT = CioCT

# Discrete time reference signal matrices
ArDT = np.matrix([[.99975, .0097951],[-.048526, .94152]])
BrDT = np.matrix([[.000024506],[.048526]])
'''

'''
# Converting matrices to discrete using Runge-Kutta
ADT = np.identity(2) + ACT*h + np.power(ACT,2)*np.power(h,2)/2 + np.power(ACT,3)*np.power(h,3)/6 + np.power(ACT,4)*np.power(h,4)/24
BDT = BCT*h + ACT*BCT*np.power(h,2)/2 + 1/6*np.power(ACT,2)*BCT*np.power(h,3) + 1/24*np.power(ACT,3)*BCT*np.power(h,4) 
CDT = CCT
DDT = DCT
'''

Phi_Phi, Phi_F, Phi_R, Ae, Be, Ce, k, u1, y1 = mpcfunc.mpcgain(ADT, BDT, CDT, Nc, Np, ref, w_f)
#Phi_Phi, Phi_F, Phi_R, Ae, Be, Ce, k, u1, y1 = mpcfunc.mpcgain(AioDT, BioDT, CioDT, Nc, Np, ref, w_f)


# Packing matrices into a state-space structure
#sys1 = signal.StateSpace(ACT, BCT, CCT, DCT)
#sys1 = sys1.to_discrete(.02)
#t1, y1, x = signal.dlsim(sys1, u=np.array([5]*100))


#plt.plot(t1, y1)

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.step(k[0,:], y1[0,:])
ax2.step(k[0,:], u1[0,:])
ax1.set_ylabel('Output')
ax2.set_xlabel('Sampling Instant')
ax2.set_ylabel('Control')
#plt.plot(k[0,:], y1[0,:], label='Output')
#plt.plot(k[0,:], u1[0,:], label='Control')
#plt.ylabel()
#plt.xlabel('Sampling Instant')
#plt.legend()
#plt.show()

