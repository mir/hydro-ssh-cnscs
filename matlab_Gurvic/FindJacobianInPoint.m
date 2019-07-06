function [jac] = FindJacobianInPoint(acostheta)
%% Modeling parameters for generator
global omega_s r T_r E_r T_rd T_rq T_J U_nom U theta_net A ...
	   x_d x_q x_ad x_r x_rd x_rq x_aq gamma;
gamma     (gamma     == -1) = 1;
U_nom     (U_nom     == -1) = 15.75 * 10^3; 
theta_net (theta_net == -1) = 0;
omega_s   (omega_s   == -1) = 142.8 * 2 * pi / 60;
T_J       (T_J       == -1) = 25.50 * 10^6; 
E_r       (E_r       == -1) = 530;
r         (r         == -1) = 0.0034;
T_r       (T_r       == -1) = 8.21;
x_d       (x_d       == -1) = 1.58;
x_q       (x_q       == -1) = 0.97;

U = U_nom * gamma; 
ddT_d = 0.143;
ddT_q = 0.243;
dx_d  = 0.43;
ddx_d = 0.3;
ddx_q = 0.31;
x_s   = 0.184;
vartheta_d = 1;
vartheta_q = 1;
x_ad =  x_d - x_s;
x_aq =  x_q - x_s;
x_r =  x_ad^2/(x_d - dx_d);
x_sr = x_r - x_ad;
x_srd = 1 /( 1 / (ddx_d-x_s) - 1/x_ad - 1/x_sr);
x_srq = 1 /( 1 / (ddx_q-x_s) - 1/x_aq);
x_rd =  x_ad + x_srd;
x_rq =  x_aq + x_srq;
r_rd = (x_rd*x_d-x_ad^2) * ddx_d / (omega_s * x_d * dx_d * ddT_d);
r_rq = (x_rq*x_q-x_aq^2) / (omega_s * x_q * ddT_q);
T_rd = x_rd / (omega_s * r_rd);
T_rq = x_rq / (omega_s * r_rq);

a11 = (x_s + vartheta_d * x_ad);
a13 = vartheta_d;
a15 = vartheta_d;
a22 = (x_s + vartheta_q * x_aq);
a24 = - vartheta_q;
a31 = vartheta_d * x_ad^2/x_r;
a33 = vartheta_d * x_ad/x_r + x_sr/x_r;
a35 = vartheta_d * x_ad/x_r;
a41 = vartheta_d * x_ad^2/x_rd;
a43 = vartheta_d * x_ad/x_rd;
a45 = vartheta_d * x_ad/x_rd + x_srd/x_rd;
a52 = vartheta_q * x_aq^2/x_rq;
a54 = - vartheta_q * x_aq/x_rq - x_srq/x_rq;

A = [ ...
    a11,   0, a13,   0, a15; ...
      0, a22,   0, a24,   0; ...
    a31,   0, a33,   0, a35; ...
    a41,   0, a43,   0, a45; ...
      0, a52,   0, a54,   0  ...    
];

%% Modeling parameters for turbune 
global S l ro p_u p_l k C Kt g Q_max;
Q_max (Q_max == -1) = 358;
S     (S     == -1) = pi / 4 * 7.5^2;
l     (l     == -1) = 196;
k     (k     == -1) = 0.43;
C     (C     == -1) = 0.1;

ro  = 0.98 * 10^3;
g   = 9.81;
p_u = ro*l*g;
p_l = 0.01 * p_u;
Kt  = S/l/ro;

%% Modeling parameters for turbune control
global CK T_0 mu_min mu_max mu_nom use_control senser_sense;

CK           (CK           == -1) = 10;
senser_sense (senser_sense == -1) = 2*pi*30*10^-3;
use_control  (use_control  == -1) = 1;
T_0          (T_0          == -1) = 10;

ttheta = acos(0.9);
mu_nom = findMu0(ttheta);
mu_max = Q_max/C/sqrt(p_u - p_l);

if (mu_nom > mu_max)
    mu_nom = mu_max;
end;

mu_min = mu_nom - mu_max;
mu_max = mu_nom - 0.05 * mu_max;

%% FIND EQUILIBRIUM POINT   
 
theta = acos(acostheta);

mm = 0;
s = 0;
Q = C*mu_nom*sqrt(p_u-p_l);

i_d = x_q/(r^2+x_d*x_q)*(r/x_q*U * sin(theta) + U * cos(theta) - E_r);
i_q =  r/(r^2+x_d*x_q)*(x_d/r*U * sin(theta) + U * cos(theta) + E_r);

psi_d  = x_d * i_d + E_r;
psi_q  = x_q * i_q;
psi_r  = x_ad^2/x_r*i_d + E_r;
psi_rd = x_ad^2/x_rd*i_d + x_ad/x_rd*E_r;
psi_rq = x_aq^2/x_rd*i_q;

jac = JacobianFunction(psi_d, psi_q, psi_r, psi_rd, psi_rq, s, theta, Q, mm);

end

function [mu_0] = findMu0(theta) 
global U E_r r x_q x_d;

u_d = @(t) - U * sin(t);
u_q = @(t)   U * cos(t);

i_d = @(t) - x_q/(r^2 + x_d*x_q) * (r / x_q * u_d(t) - u_q(t) + E_r);
i_q = @(t) - r  /(r^2 + x_d*x_q) * (x_d / r * u_d(t) + u_q(t) - E_r);

psi_d = @(t) x_d * i_d(t) + E_r;
psi_q = @(t) x_q * i_q(t);

phi = @(t) psi_d(t) * i_q(t) - psi_q(t) * i_d(t);

global k C p_u p_l omega_s;

koef = k * C * (p_u - p_l)^(3/2) / omega_s^2;

mu_0 = phi(theta) / koef;

end