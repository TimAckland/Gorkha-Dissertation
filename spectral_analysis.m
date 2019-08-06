%take one of the Rcs from SSA
RC = RC9;


hours = 1:length(RC);

%pad data - because FFT will not handle NaNs
bad_vals = isnan(RC);
good_vals = ~isnan(RC);
% RC(bad_vals) = 0.0;  %pad with zero
RC(bad_vals) = mean(RC(good_vals));  %pad with series mean


figure(5)
subplot(3,1,1)
plot(hours,RC);
title('RC')
ylabel('Variance')
xlabel('Hours')


n=length(RC);

% compute FFT using built-in function
Y = fft(RC);
Y(1)=[]; %remove first component (sum of data)

% and compute power versus frequency plot (periodogram)
power = abs(Y(1:floor(n/2))).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;  % make a frequency axis

subplot(3,1,2)
plot(freq,power)
xlabel('Cycles/hour')
ylabel('Power')

subplot(3,1,3)
% plot period instead of frequency
period=1./freq;  % period axis as inverse of frequency
plot(period,power);
axis([0 4000 0 max(power)*1.25]);
ylabel('Power');
xlabel('Period (Hours)');


% optional extra: locate and mark the peak frequency (see if you can 
% figure this out yourself!)
hold on; % this enables overplotting with deleting existing figure content
index=find(power==max(power)); % note == rather than = to test equivalence
mainPeriodStr=num2str(period(index));
plot(period(index),power(index),'r.', 'MarkerSize',25);
text(period(index)+2,power(index),['Period = ',mainPeriodStr]);
hold off; %turns off overplotting



