function [] = drawResult( t, z, color )
global theta_net;

theta  = z(:,1); s = z(:,2); mu = z(:,9);

size = 30;
lineWidth = 2;
thetaLabelText = '\theta_\Delta';
sLabelText     = 's';
muLabelText    = '\mu_\Delta';
PLabelText     = 'P_\Delta';

figure(1);
subplot(4, 1, 1), 
  hold on;
  plot(t, theta, '-', ...
      'Color', color, ...
      'LineWidth', lineWidth);

  xlabel(thetaLabelText, ...
      'FontSize', size);
  grid on;

subplot(4, 1, 2), 
  hold on;
  plot(t, s, '-', ... 
      'Color', color, ...
      'LineWidth', lineWidth);

  xlabel(sLabelText, ...
      'FontSize', size);
  grid on;

subplot(4, 1, 3), 
  hold on;
  plot(t, mu, '-', ...
      'Color', color, ...
      'LineWidth', lineWidth);

  xlabel(muLabelText, ...
      'FontSize', size);    
  grid on;

subplot(4, 1, 4), 
  hold on;
  plot(t, (findPower(theta + theta_net) ...
      - findPower(acos(0.9))) ./ 10^6, '-', ... 
      'Color', color, ...
      'LineWidth', lineWidth);

  xlabel(PLabelText, ...
      'FontSize', size);    
  grid on;

figure(2);
  hold on;
  plot3(s, theta, mu, '-', ... 
      'Color', color, ...
      'LineWidth', lineWidth);    
 
  xlabel(sLabelText, ...
      'FontSize', size);  
  ylabel(thetaLabelText, ...
      'FontSize', size);  
  zlabel(muLabelText, ...
      'FontSize', size);  
  
  grid on;
end

function [P] = findPower(theta) 
global E_r r x_d x_q U omega_s;
u_d = @(tt, uu) - uu .* sin(tt);
u_q = @(tt, uu)   uu .* cos(tt);

i_d = @(tt, uu) ...
    - x_q./(r.^2 + x_q.*x_d) ...
    .*(r ./ x_q.*u_d(tt, uu) ...
     - u_q(tt, uu) + E_r);
i_q = @(tt, uu) ...
    - r  ./(r.^2 + x_q.*x_d) ...
    .*(x_d ./ r.*u_d(tt, uu) ...
    + u_q(tt, uu) - E_r);

psi_d = @(tt, uu) x_d ...
            .* i_d(tt, uu) + E_r;
psi_q = @(tt, uu) x_q ...
            .* i_q(tt, uu);

p = @(tt, uu) omega_s ...
    .* (i_q(tt, uu).*psi_d(tt, uu) ...
    - i_d(tt, uu).*psi_q(tt, uu));

P = p(theta, U);
end