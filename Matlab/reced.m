clear all; clc; close all;

Bp = [.0012171; .24016];
Ap = [.99878 .0096065;-.24016 .92193];
Cp = [1 0];
% Ap = 0.6;
% Bp = 1;
% Cp = 1;
Dp = 0;
Nc = 4;
Np = 10;
w_f = 20;
ref = 5;

[Phi_Phi, Phi_F, Phi_R, A_e, B_e, C_e, Phi, F]=mpcgain(Ap, Bp, Cp, Nc, Np);

[n, n_in]=size(B_e);
xm=[0; 0];
Xf = zeros(n, 1);
N_sim=400;
r=ones(N_sim,1)*ref;
u=0;
y=0;

for kk=1:N_sim
    DeltaU=inv(Phi_Phi+w_f*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
    %DeltaU=inv(Phi_Phi+w_f*eye(Nc,Nc))*Phi'*(r(kk)-F*Xf);
    deltau=DeltaU(1,1);
    u=u+deltau;
    u1(kk)=u;
    y1(kk)=y;
    xm_old=xm;
    xm=Ap*xm+Bp*u;
    y=Cp*xm;
    Xf=[xm-xm_old;y];
end

k=0:(N_sim-1);
figure
subplot(211)
plot(k,y1)
xlabel('Sampling Instant')
ylabel('Output')
%legend('Output')
subplot(212)
plot(k,u1)
xlabel('Sampling Instant')
ylabel('Control')
%legend('Control')


