% Jason G. Pickel
% March 6 2017
% Determining the linearizing PHX dynamics about steady-state

clc; clear all; close all;

format compact
format shortg

%%%%%%%% Parameters
hh = sym('hh','real');
ch = sym('ch','real');
A = sym('A','real');
Vh = sym('Vh','real');
rhoh = sym('rhoh','real');
mdoth = sym('mdoth','real');
Thi   = sym('Thi','real');
hc = sym('hc','real');
cc = sym('cc','real');
Vc = sym('Vc','real');
rhoc = sym('rhoc','real');
u = sym('u','real');
Tci  = sym('Tci','real');
wt   = sym('wt','real');
cw = sym('cw','real');
A = sym('A','real');
Vw = sym('Vw','real');
wp = sym('wp','real');
rhow = sym('rhow','real');
x1bar = sym('x1bar','real');
x2bar = sym('x2bar','real');
x3bar = sym('x3bar','real');
x4bar = sym('x4bar','real');
mdotcbar = sym('mdotcbar','real');
b1    = sym('b1','real'); %variation with respect to SS
b2    = sym('b2','real');
b3    = sym('b3','real');
b4    = sym('b4','real');

%%%%% State coordinate translation
x1 = b1 + x1bar; %hot leg temp out
x2 = b2 + x2bar; %cold leg temp out
x3 = b3 + x3bar; %wall temp
x4 = b4 + x4bar; %cold leg mass flow rate out of pump
mdotc = u + mdotcbar;

%%%%% State Dynamics
f1 = -(mdoth*ch + hh*A)/(rhoh*Vh*ch)*x1 + hh*A/(rhoh*Vh*ch)*x3...
    + mdoth*ch/(rhoh*Vh*ch)*Thi;
f2 = -(x4*cc + hc*A)/(rhoc*Vc*cc)*x2 + hc*A/(rhoc*Vc*cc)*x3...
    + x4*cc/(rhoc*Vc*cc)*Tci*(1+wt);
f3 = hh*A/(rhow*Vw*cw)*x1 + hc*A/(rhow*Vw*cw)*x2...
    - (hh*A + hc*A)/(rhow*Vw*cw)*x3;
f4 = -wp*x4 + wp*mdotc;

%%%%% Determining the Linearizing Dynamics
f1_tay = taylor(f1,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',2);
f1_bias = taylor(f1,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',1);
f1_lin  = f1_tay - f1_bias;
A1      = gradient(f1_lin,[b1,b2,b3,b4])';
B1w      = gradient(f1_lin,wt);
B1u      = gradient(f1_lin,u);

f2_tay = taylor(f2,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',2);
f2_bias = taylor(f2,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',1);
f2_lin  = f2_tay - f2_bias;
A2      = gradient(f2_lin,[b1,b2,b3,b4])';
B2w      = gradient(f2_lin,wt);
B2u      = gradient(f2_lin,u);

f3_tay = taylor(f3,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',2);
f3_bias = taylor(f3,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',1);
f3_lin  = f3_tay - f3_bias;
A3      = gradient(f3_lin,[b1,b2,b3,b4])';
B3w     = gradient(f3_lin,wt);
B3u      = gradient(f3_lin,u);

f4_tay = taylor(f4,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',2);
f4_bias = taylor(f4,[b1,b2,b3,b4,u],[0,0,0,0,0],'Order',1);
f4_lin  = f4_tay - f4_bias;
A4      = gradient(f4_lin,[b1,b2,b3,b4])';
B4w      = gradient(f4_lin,wt);
B4u      = gradient(f4_lin,u);

%%%%% Linearized Matrices
% State matrix
disp('State matrix')
A=[A1;A2;A3;A4]

% Input matrix noise inputs
disp('Input matrix noise')
Bw = [B1w;B2w;B3w;B4w]

% Input matrix control input
disp('Input matrix control')
Bu = [B1u;B2u;B3u;B4u]













