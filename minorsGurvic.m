function [ res ] = minorsGurvic()
jac = JacobianFunction(0);

syms lambda;
d = - det(jac - lambda * eye(9));
a=[];
for i=0:9
 a = [subs(diff(d, i) / factorial(i), ...
     lambda, 0), a];
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

