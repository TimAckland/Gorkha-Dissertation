function [xmean,xstd,xk,ro,lam]=vssanan(ts,m,f,raw_variance)
%function [xmean,xstd,xk,ro,lam]=vssanan(ts,m,f, raw_variance) 
% 
% purpose: given time series ts and window size m (m*dt), compute 
%          principal components xk, eigenvectors ro, and eigenvalues lam. 
%          see Vautard and Ghil 1989 for notation 
%          this version ignores NaNs (i.e. handles missing data in series) 
% 
% inputs:  ts=time series with constant dt 
%          m=window size, in time steps 
%          f=fraction (0<f<=1) of good data points for determining PCs 
%          raw_variance=variance of series from earlier SSA passes
% outputs: xmean=mean of ts 
%          xstd=standard deviation of ts 
%          xk=principal components zero mean time series, row k = pc k 
%          ro=eigenvectors ro(j,k), k=corresponding PC 
%          lam=eigenvalue vector, sorted from greatest to smallest 
% 
% by David Schoellhamer, 5/27/93, 12/99, 1/2000 
% modified by JF July 2007

disp('Performing SSA with allowance for missing data (Schoellhamer method)')

% normalize time series x 
igood=find(~isnan(ts)); 
xmean=mean(ts(igood)); 
xstd=std(ts(igood)); 
x=(ts-xmean)/xstd; 
n=length(x);

%JF change:
%if original series variance from earlier passes specified, derive
%correction factor to be applied to eigenvalues
if nargin == 4
 correction =   (xstd*xstd)/raw_variance;
else
 correction = 1;
end


% compute autocorrelation function eqn b.2 
disp(' ')
disp('Computing ACF') 
% only need the first m values 
c = zeros(m,1);%JF change
for j1=1:m 
%  summ=0; 
  j=j1-1; 
%  for i=1:n-j 
%    summ=summ+x(i)*x(i+j); 
%  end 
%  c(j1)=summ/(n-j-1); 
% 
% modify for nan 
  prod=x(1:n-j).*x(j+1:n); 
  igood=find(~isnan(prod)); 
  c(j1)=sum(prod(igood))/(length(igood)-1); 
end 

% form covariance matrix a, divide by m (eqn b.3) 
disp('Forming covariance matrix a')  
a=toeplitz(c)/m; 
%JF change
disp([' Total variance (should == 1 for normalised data): ' num2str(trace(a))])

% determine eigenvalues and eigenvectors of a (solve b.3) 
disp('Determining eigenvalues and eigenvectors'); 
[z,eval]=eig(a); 
lam=diag(eval); 
%JF change
lam = lam .* correction; %Apply any correction to variance due to prior SSA


% sort by descending eigenvalues 
disp('Sorting by descending eigenvalues'); 
disp(' ')
[lam,ilam]=sort(lam); 
lam=flipud(lam); 
ilam=flipud(ilam); 
summ=0; 
kmax=0;
%JF changes here
order = 1:length(lam); order = order';
error_interval(order) = 1.96 .* lam(order) .* sqrt(m/(2*n));error_interval = error_interval';
figure(10)
% semilogy(lam(1:min([m 30])),'o')
% plot(lam(1:min([m 30])),'o')
max_plot_lag = 30; %usual default max no of lags to plot
% max_plot_lag = 60;
errorbarlog(order(1:min([m max_plot_lag])),lam(1:min([m max_plot_lag])),error_interval(1:min([m max_plot_lag])),'o');
nkmax=0; 
while nkmax==0 
  nkmax=input('Number of pcs (0 for keyboard control) = '); 
  if nkmax==0 
    disp('Enter return command to resume') 
    keyboard 
  end 
end 
while kmax<nkmax 
  kmax=kmax+1; 
  ro(:,kmax)=z(:,ilam(kmax)); 
  summ=summ+lam(kmax); 
  disp(sprintf('Eigenvalue of pc %d is %g, sum of variance is %g',... 
               kmax,lam(kmax),summ)) 
end 

% determine principal component time series (eqn b.6) 
% 
xk=nan*ones(kmax,n-m+1); 
disp(' ');
for k=1:kmax 
%   disp(sprintf('Determining principal component time series %d',k))
  fprintf('Determining principal component time series %d\n',k)
  for i=0:n-m 
%   modify for nan 
    prod=x(i+1:i+m).*ro(:,k); 
    igood=find(~isnan(prod)); 
    ngood=length(igood); 
%   must have at least m*f good points 
    if ngood>=m*f 
      xk(k,i+1)=sum(prod(igood))*m/ngood; 
    end 
  end 
end 

