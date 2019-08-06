function [Rk]=reconstructssa(N,M,xmean,xstd,xk,ro,kset) 
%function [Rk]=reconstructssa(N,M,xmean,xstd,xk,ro) 
%         [Rk]=reconstructssa(N,M,xmean,xstd,xk,ro,kset)
%
% purpose:    given SSA results from SSANAN.M, reconstruct 
%             time-series using kmax principal components drawn from xk 
%             see Vautard and Ghil 1992 for algorithm 
%             This function by JF 5/2007 
% 
% inputs:     N=length of original time series 
%             M=window size, in time steps
%             xmean=mean of original time-series
%             xstd=standard deviation of original time-series
%             xk=principal components zero mean time series, row k = pc k 
%             ro=eigenvectors ro(j,k), k=corresponding pc
%(optionally) kset = vector of PCs to use
%
% outputs: Rk = reconstructed series (from k PCs) of length N 
%

if nargin == 7 %we have passed the number (or a set) of RC to draw from xk
  kmax = length(kset);
  %Pull out required rows of xk and columns of ro
  xk = xk(kset,:);
  ro = ro(:,kset);
end
if kmax > min(size(xk))
 disp(['Warning: kmax > min(size(xk)) : reset to ' numstr(min(size(xk)))]);
end 

%Reconstruct series from k PCs according to Vautard et al (1992)
reciprocal_M = 1 / M;
Rk = zeros(N,1);
%Central portion of series
for t = M:N-M+1
 product = 0;
 for k = 1:kmax
  for j = 1:M
   product = product + xk(k,(t-j+1)) * ro(j,k);
  end
 end
 Rk(t) = reciprocal_M * product;
end
%Start of series
for t = 1:M-1
 product = 0;
 for k = 1:kmax
  for j = 1:t
   product = product + xk(k,(t-j+1)) * ro(j,k);
  end
 end
 Rk(t) = (1/t) * product;
end
%End of series
for t = N-M+2:N
 product = 0;
 for k = 1:kmax
  for j = t-N+M:M
   product = product + xk(k,(t-j+1)) * ro(j,k);
  end
 end
 Rk(t) = (1/(N-t+1)) * product;
end

% Rk = Rk * xstd + xmean; %remove normalisation applied in VSSANAN
