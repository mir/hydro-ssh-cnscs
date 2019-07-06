function [ ] = setModelingParameters(gK)
%% simulateModel sets paremeters 
%  for system and starts modelling
%  
%  parameter g     is koefficient 
%  parameter time  is modelling time
%  parameter z0    is initial data

global omega_s r T_r E_r T_rd T_rq ...
       T_J U_nom U theta_net A x_d ...
       x_q x_ad x_r x_rd x_aq gamma ...
       S l ro p_u p_l k C Kt g Q_max ...
       CK T_0 mu_min mu_max mu_nom ...
       use_control senser_sense; 

%%  Modeling parameters for generator
    gamma = gK; U_nom = 15.75 * 10^3; 
    theta_net = acos(0.9); 
    omega_s = 142.8 * 2 * pi / 60; 
    T_J = 25.50 * 10^6; E_r = 530;
    r = 0.0034; T_r = 8.21; x_d = 1.58; 
    x_q = 0.97; ddT_d = 0.143; ddT_q = 0.243;
    dx_d = 0.43; ddx_d = 0.3; ddx_q = 0.31; 
    x_s = 0.184; U = U_nom * gamma;
    x_ad =  x_d - x_s; x_aq = x_q - x_s;
    x_r =  x_ad^2/(x_d - dx_d); 
    x_sr = x_r - x_ad;
    x_srd = 1 ...
         /( 1 / (ddx_d-x_s) - 1/x_ad - 1/x_sr);
    x_srq = 1 ...
         /( 1 / (ddx_q-x_s) - 1/x_aq);
    x_rd =  x_ad + x_srd; x_rq =  x_aq + x_srq;
    r_rd = (x_rd*x_d-x_ad^2) * ddx_d ...
           / (omega_s * x_d * dx_d * ddT_d);
    r_rq = (x_rq*x_q-x_aq^2) ...
           / (omega_s * x_q * ddT_q);
    T_rd = x_rd / (omega_s * r_rd); 
    T_rq = x_rq / (omega_s * r_rq);

    A = [ x_d,0,1,0,1;...
        0,x_q,0,-1,0;...
        x_ad^2/x_r,0,1,0,x_ad/x_r;...
        x_ad^2/x_rd,0,x_ad/x_rd,0,1;...
        0,x_aq^2/x_rq,0,-1,0];

%% Modeling parameters for turbune 
    Q_max = 358; S = pi / 4 * 7.5^2; 
    l = 196; k = 40; C = 0.1; 
    ro = 0.98 * 10^3; g = 9.81; 
    p_u = ro*l*g; p_l = 0.01 * p_u;
    Kt = S/l/ro;

%% Modeling parameters for turbune control
    CK = 10; senser_sense = 2*pi*30*10^-3; 
    use_control = 1; T_0 = 10;
    
    u_d = @(t) - U * sin(t); 
    u_q = @(t)   U * cos(t);
    i_d = @(t) - x_q/(r^2 + x_d*x_q) ...
        * (r / x_q * u_d(t) - u_q(t) + E_r);
    i_q = @(t) - r  /(r^2 + x_d*x_q) ...
        * (x_d / r * u_d(t) + u_q(t) - E_r);
    psi_d = @(t) x_d * i_d(t) + E_r; 
    psi_q = @(t) x_q * i_q(t);
    phi = @(t) psi_d(t) ...
     * i_q(t) - psi_q(t) * i_d(t);

    mu_nom = phi(acos(0.9)) ...
            / (k * C * (p_u - p_l)^(3/2) ...
                / omega_s^2);
    mu_max = Q_max/C/sqrt(p_u - p_l);


    if (mu_nom > mu_max)
        mu_nom = mu_max;
    end;

    mu_min = mu_nom - mu_max;
    mu_max = mu_nom - 0.05 * mu_max;
end