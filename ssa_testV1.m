%Experiment with analysis of groundwater series SSA

%Test dataset is hourly series of groundwater (via Willy Burgess)
%Site: Laksmipur PZ1
%X is time (hours, arbitrary origin)
%Y is groundwater head

clear
close all

%edit these values as required
M = 250; %Window size (in N steps of dt)
f = 1; %Fraction of good data (range 0 to 1; assume 1 if series complete);


%specify input file
[filename,pathname] = ...
  uigetfile('*.csv','Specify input file (2 column CSV format)');
  infile = fullfile(pathname,filename);

disp(' ')
disp(' ')
disp('SSA analysis')
disp(['Data =  ' filename])


%load data file
fid = fopen(filename,'r');

%skip the header
header_line = fgetl(fid);

%pre-dimension arrays larger than we need 
X = ones(100000,1);
Y = X;

i = 0;
while ~feof(fid)
 i = i + 1;
 line = fgetl(fid);
 comma = strfind(line,',');
 X(i) = str2double(line(1:comma-1));
 Y(i) = str2double(line(comma+1:end));
end

fclose(fid);

N = i;

%trim the X and Y vectors to the data length
X = X(1:N);
Y = Y(1:N);

%series to be analysed is Y, with unit step (= 1 hour)
series_variance = std(Y)^2;

%Initial QC
figure(1)
subplot(2,1,1)
 plot(X,Y)

%look for spikes or dropouts
minY = 1;  %edit as appropriate
maxY = 130000;

dropouts = find(Y < minY);
spikes = find(Y > maxY);
if ~isempty(dropouts)
 Y(dropouts) = NaN;
end
if ~isempty(spikes)
 Y(spikes) = NaN;
end

subplot(2,1,2)
 plot(X,Y)



%Perform SSA
%Call vssanan.m (slightly modified from Schoellhamer (2001) code)
%This will prompt user for number of PCs to evaluate
[xmean,xstd,xk,ro,lam] = vssanan(Y,M,f);

%Reconstruct series from major PCs according to Vautard et al (1992)
%Call reconstructssa.m (French, 2007)
s = size(ro); N_RC = s(2);
% N_RC = 5; %just use first 5 RCs
Rk = reconstructssa(N,M,xmean,xstd,xk,ro,1:N_RC); %based on all PCs
Rk = Rk * xstd + xmean; %remove normalisation applied in VSSANAN

RC1 = reconstructssa(N,M,xmean,xstd,xk,ro,[1]); 
RC2 = reconstructssa(N,M,xmean,xstd,xk,ro,[2]); 
RC3 = reconstructssa(N,M,xmean,xstd,xk,ro,[3]); 
RC4 = reconstructssa(N,M,xmean,xstd,xk,ro,[4]); 
RC5 = reconstructssa(N,M,xmean,xstd,xk,ro,[5]); 
RC6 = reconstructssa(N,M,xmean,xstd,xk,ro,[6]); 
RC7 = reconstructssa(N,M,xmean,xstd,xk,ro,[7]); 
RC8 = reconstructssa(N,M,xmean,xstd,xk,ro,[8]); 
RC9 = reconstructssa(N,M,xmean,xstd,xk,ro,[9]); 
% RC10 = reconstructssa(N,M,xmean,xstd,xk,ro,[10]); 
% RC11 = reconstructssa(N,M,xmean,xstd,xk,ro,[11]);

Num = 10
figure(3)
subplot(Num,1,1)
plot(Rk,'b') %reconstructed series based on 'n' PCs
title('Reconstructed');
%following plots are for specific RCs - edit as required
subplot(Num,1,2)
plot(RC1)
title('RC1');
subplot(Num,1,3)
plot(RC2)
title('RC2');
subplot(Num,1,4)
plot(RC3)
title('RC3');
subplot(Num,1,5)
plot(RC4)
title('RC4');
subplot(Num,1,6)
plot(RC5)
title('RC5');
subplot(Num,1,7)
plot(RC6)
title('RC6');
subplot(Num,1,8)
plot(RC7)
title('RC7');
subplot(Num,1,9)
plot(RC8)
title('RC8');
subplot(Num,1,10)
plot(RC9)
title('RC9');
% subplot(Num,1,11)
% plot(RC10)
% title('RC10');
% subplot(Num,1,12)
% plot(RC11)
% title('RC11');


