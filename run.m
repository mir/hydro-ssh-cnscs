z = [0 1 0 0 0 0 0 0 0];

setModelingParameters(1);

time = [0, 50];
[t, z] = ode45(@model, time, z(end,:));

t_span = 500;
for i=1:1
time = [time(end), time(end) + t_span];
[t, z] = ode45(@model, time, z(end,:));
drawResult(t - 500, z, [0.1*i, 0, 1 - 0.1*i]);
end

setModelingParameters(0.89);

for i=2:10
time = [time(end), time(end) + t_span];
[t, z] = ode45(@model, time, z(end,:));
drawResult(t - 500, z, [0.1*i, 0, 1 - 0.1*i]);
end