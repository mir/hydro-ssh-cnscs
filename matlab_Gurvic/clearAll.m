function [ ] = clearAll()
  
clearGlobalProperties();
    
end

function [] = clearGlobalProperties()
global  omega_s gamma r T_r E_r T_J U_nom theta_net x_d x_q S l k C CK T_0 ...
        mu_min mu_max mu_nom use_control senser_sense acos_val Q_max;

acos_val = -1; Q_max = -1; gamma = -1; omega_s = -1; r = -1; T_r = -1; 
E_r = -1; T_J = -1; U_nom = -1; theta_net = -1; x_d = -1; x_q = -1; S = -1; 
l = -1; k = -1; C = -1; CK = -1; T_0 = -1; mu_min = -1; mu_max = -1; 
mu_nom = -1; use_control = -1; senser_sense = -1; 

end

