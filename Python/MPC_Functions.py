
import numpy as np

def mpcgain(Ap, Bp, Bw, Cp, Nc, Np, ref, w_f, N_sim):
    # Compute gain matrices
    if len(Cp.shape) == 1:
        n1 = len(Cp)
        m1 = 1
    else:
        m1 = Cp.shape[0]
        n1 = Cp.shape[1]
    n_in = Bp.shape[1]
    Ae = np.identity(n1+m1)
    Ae[0:n1,0:n1] = Ap
    Ae[n1:n1+m1, 0:n1] = np.dot(Cp,Ap)
    Be = np.zeros((n1+m1,n_in))
    Be[0:n1,:] = Bp
    Be[n1:n1+m1,:] = np.dot(Cp, Bp)
    Ce = np.zeros((m1, n1+m1))
    Ce[:,n1:n1+m1] = np.identity(m1)
    n = n1+m1
    h = np.zeros((Np, Ce.shape[1]))
    F = np.zeros((Np, np.dot(Ce,Ae).shape[1]))
    h[0,:] = Ce
    F[0,:] = np.dot(Ce, Ae)
    for kk in range(1,Np):
        h[kk,:] = np.dot(h[kk-1,:], Ae)
        F[kk,:] = np.dot(F[kk-1,:], Ae)
    v = np.dot(h, Be)
    Phi = v
     
    for i in range(1, Nc):
        w = v[0:Np-i,0][np.newaxis]
        w = np.transpose(w)
        append2Phi = np.vstack((np.zeros((i,1)),w))
        Phi = np.hstack((Phi, append2Phi))

    BarRs=np.ones((Np, 1))
    Phi_Phi = np.dot(np.transpose(Phi), Phi)
    Phi_F = np.dot(np.transpose(Phi), F)
    Phi_R = np.dot(np.transpose(Phi),BarRs)
    
    #Carry out receding horizon control
    n2 = Be.shape[0]
    xm = np.zeros((n1,1))
    Xf = np.zeros((n2,1))
    w = np.random.normal(0,.1,N_sim)
    r = ref*np.ones((N_sim, 1))
    u = 0
    y = 0
    u1 = np.zeros((1, N_sim))
    y1 = np.zeros((1, N_sim))
    
    for kk in range(1, N_sim):
        DeltaU = np.dot(np.linalg.inv(Phi_Phi+w_f*np.identity(Nc)),Phi_R*r[kk]-np.dot(Phi_F,Xf))
        deltau = DeltaU[0,0]
        u = u + deltau
        u1[0,kk] = u
        y1[0,kk] = y
        xm_old = xm
        #print(w[kk])
        #print(xm)
        xm = np.dot(Ap,xm) + np.dot(Bp,u) #+ np.dot(Bw, w[kk])
        y = np.dot(Cp,xm)
        #print(y)
        Xf = np.vstack((xm-xm_old,y))
        
    k = np.linspace(0,N_sim-1, N_sim)
    k = k.reshape((1,N_sim))

    return Phi_Phi, Phi_F, Phi_R, Ae, Be, Ce, k, u1, y1
    
def QPHild(E, F, M, gamma):
    # Determine which constraints are active and which are inactive
    n1 = M.shape[0]
    m1 = M.shape[1]
    x = np.dot(-np.linalg.inv(E),F)
    kk=0
    
    '''
    for i in range(1,n1):
        if M[i,:]*x > gamma[i]:
            kk += 1
        else:
            kk = kk+0
    '''
    H = np.dot(M,np.dot(np.linalg.inv(E),np.transpose(M)))
    K = np.dot(M,np.dot(np.linalg.inv(E),F)) + gamma
    #print(H)
    #print(K)
    n = K.shape[0]
    m = K.shape[1]
    x_ini = np.zeros((n, m))
    lam = x_ini
    al = 10;
            
    for km in range(1,400):
        lambda_p = lam
        
        for i in range(0,n):
            #print(i)
            #print(lambda_p)
            w = np.dot(H[i,:],lam) - np.dot(H[i,i],lam[i,0])
            w = w + K[i,0]
            la = -w/H[i,i]
            #print(la)
            #print(w)
            lam[i,0] = max(0,la)
            #print(lam[i,0])
        al = np.dot(np.transpose(lam-lambda_p),lam-lambda_p)
        #print(lambda_p)
        #print(al)
        #print(lam)
        if al < 10e-8:
            break
    #print(km)
    x = x - np.dot(np.linalg.inv(E),np.dot(np.transpose(M),lam))

    return x



        