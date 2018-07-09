function [adj_temp,seg_temp,fs,adj_time] = shallow_model(Time,temp,zz,K,Gz)
% This matlab function can calculate a temperature record at a specific seafloor 
% sediment depth based on seafloor/ocean bottom temperature record.
%
% This m-file uses FFT to deconstruct the temperature record and the equation from
% Hamamato et al, 2005 to reconstruct the temperature profile.
%
% Inputs:
%	
%	Time: An array of the time (units of days) 
%       or the serial matlab time format (see datenum for more information)
%
%	temp: An array of seafloor temperatures (degree C). Matrix should be the same size as Time
%	zz: sediment depth (meters) that the temperature will be propagated to
%	K: Sediment thermal diffusivity (m^{2} s^{-1})
%	Gz: Sediment geothermal gradient (degree C m^{-1})
%
% Outputs: 
% 
%	adj_temp: Estimated temperatures at depth zz
%	seg_temp: Interpolated seafloor temperature record to match adj_temp size
%	fs: sample frequency (Hz)
%	adj_time: Interpolated matrix of time to match adj_temp size
%
%   Example of usage: [adj_temp,seg_temp,fs,adj_time] = shallow_model(Time,temp,0.75,4e-7,0.105);
% 
% Marie S. Salmi - 2018

if nargin<1
	help shallow_model
	adj_temp= [];
    	seg_temp = [];
    	fs = [];
    	adj_time = [];
    	return
end

if length(Time) ~= length(temp)
	disp(' The time and temperture matrices are not the same length')
    	adj_temp= [];
    	seg_temp = [];
    	fs = [];
	adj_time = [];
	return
end

%Sample matrix length
L = length(Time);
a=nextpow2(L);

% Reduce the data size if the input is extremely large
if L > 50000
    N=2^(a-2);
else
    N=2^(a);
end

%load and interpolate the data
min_time=min(Time);
max_time=max(Time);

tp=linspace(min_time,max_time,N)-min_time;
Temp=interp1(Time-min_time,temp,tp');

%variables used in the calculation. Includes Depth and K. 
seg = 10;                 %Used to reduce the dataset to speed up the processing
N2  = (N)/2;                                  %Sample cutoff
fs  = (1/(tp(1,2)-tp(1,1)))*(1/(60*60*24));    %sample frequency in seconds
Tz  = Temp-mean(Temp);

%Preallocate matrices for calculations
t=1;
rho = zeros(1,N2);
b = zeros(1,N2);
complex = zeros(1,N2);
z = zeros(1,N2);

%Run the FFT and separate the result
fast=fft(Tz,N);
realz=real(fast);
img=imag(fast);

%Actual calculations    
%h= waitbar(0, 'Calculating....');
while t <= N
    
    for n = 1:N2
        rho(n)=(2*pi*(t-1))*((n-1)/N);
        b(n)=(-(sqrt(((2*pi*(n-1))/N)*(fs/K)))*zz);
        complex(n)=((realz(n)*cos(rho(n)+b(n)))-(img(n)*sin(rho(n)+b(n))))*(exp(b(n)));
    end
    
    z(t)=(sum(complex))*(2/(N));
    %waitbar(t/N,h)
    t=t+seg;
end

%close(h);
Y=length(Temp);
adj_temp=z(1:seg:Y)+mean(Temp)+Gz*zz;
seg_temp=Temp(1:seg:Y);
adj_time=tp(1:seg:Y);

%Plot the results
%t=1:1:length(adj_temp);
%plot(t,adj_temp,'green')
%hold on
%plot(seg_temp,'black') 
