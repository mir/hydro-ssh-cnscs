function [ res ] = minorsGurvic(g)

% psi_d psi_q psi_r psi_rd psi_rq s theta Q mu
clearAll();

global  omega_s gamma r T_r E_r T_J U_nom theta_net x_d x_q S l k C CK ...
        T_0 mu_min mu_max mu_nom use_control senser_sense Q_max;

k = 40;
gamma = g;

acostheta = 0.9;

jac = FindJacobianInPoint(acostheta);

syms lambda;

d = - det(jac - lambda*eye(9));

a=[];

for i=0:9
    a = [subs(diff(d, i) / factorial(i), lambda, 0), a];
end

Gurvic = zeros(9);

for i=1:9
    for j=1:9
        if (i+1+2*(j-i) < 11 && i+1+2*(j-i) > 0)
            Gurvic(i,j) = a(2*(j-i)+i+1);
        end
    end
end

res = 0;

for i=1:9
    res = res + (det(Gurvic(1:i,1:i))>0);
end

end

