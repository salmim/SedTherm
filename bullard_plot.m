function [HF,r2,TG]= bullard_plot(temp,k,depth)
% This matlab function calculates the heat flow using the Bullard plot method from
% from (co-registered) temperature and thermal conductivity measurements at known
% depths into the sediment.
%
%    NOTE: This function only calculates heat flow on ONE profile at a time
% 
% Inputs-
% 	temp: Temperature measurements (Degree C)
% 	k: Thermal Conductivity (W m^{-1}K^{-1})
% 	depth: depth of the temperature/thermal conductivity measurements (meters)
%
% temp, k, and depth should be the same size
%
% Outputs-
% 	HF: Heat Flow (W m^(-2))
% 	r2: Coefficient of determination or R^2
% 	TG: Thermal Gradient (degree C meter^(-1))
%

% check that the data is in the correct format
M = size(temp);
if M(1,1) ~= 1
    temp = temp';
end

N = size(depth);
if N(1,1) ~= 1
    depth = depth';
end

L = size(k);
if L(1,1) ~= 1
    k = k';
end

if length(depth) == length(temp) && depth(1) ~= 0
	depth = [0,depth];
end

% start of calculations
del_z(1) = depth(1);

for i = 2:(length(depth))
    del_z(i) =depth(i)-depth(i-1);
end

dep = del_z(2:end);

% Calculate the thermal resistance
rest = dep(1)./k(1) + cumtrapz(dep./k);

% Determine how well the linear line fits the data
[p,~] = polyfit(temp,rest,1);
polydata = polyval(p,temp);
res = sum((rest- polydata).^2);
tot = sum((rest - mean(rest)).^2);
r2 = 1-(res/tot);           

% calculate the heat flow
HF = 1/p(1);

% calculate the thermal gradient
[q,~] = polyfit(temp,depth(2:end),1);
TG = 1/q(1);

% plot the input data and the fitted line
figure(3)
plot(temp,rest,'.-','LineWidth',2,'MarkerSize',14)
hold on
plot(temp,(temp*p(1)+p(2)),'r','LineWidth',1.5)
ylabel('Thermal Resistance (W m^{-1}K^{-1})')
xlabel('Temperature (\circC)')
grid on
axis ij
hold off
