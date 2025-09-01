%Dem0a.m
% Needs  dwtr.m and idwtr.m
clear all;close all force; clc
disp('Dem 0a: Test dwtr/idwtr -- Stand Alone Wavelet m-files')
 lw = 2.5; 
 set(0, 'DefaultAxesFontSize', 16);
 fs = 15;
 msize = 10;
 %-------------------------------------------------------------
N=1024;  %length of the signal
t=linspace(0,1,N);
% Celebrated Doppler!
data = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05));
%Symmlet 4 (8 tap filter)
hfilt = [ -0.07576571478934  -0.02963552764595  ...
          0.49761866763246   0.80373875180522  ...
          0.29785779560554  -0.09921954357694  ...
         -0.01260396726226   0.03222310060407];  
         
wavdata = dwtr(data, 5, hfilt); %WT of data by 'filt' @ 3 levels of detail

s = figure('Renderer', 'painters', 'Position', [10 10 1500 1500]);
subplot(6,1,1)
plot(data)
axis([1 N -2 2])
yline(-1.0, 'b-.', 'LineWidth', 2); yline(1.0, 'b-.', 'LineWidth', 2); grid on; %axis tight
title('Original Doppler Signal')


subplot(6,1,2)
wavdata = dwtr(data, 1, hfilt); %WT of data by 'filt' @ 3 levels of detail
plot(wavdata, '-');
hold on; 
plot(wavdata(1:512), 'k-', 'LineWidth', 2);
plot([512.5 512.5],[-2,2],'r-', 'LineWidth', 2);
yline(-1.0, 'b-.', 'LineWidth', 2); yline(1.0, 'b-.', 'LineWidth', 2); grid on; 
title('1st Decomposition (j = 1) ')
axis tight

subplot(6,1,3)
wavdata = dwtr(data, 2, hfilt); %WT of data by 'filt' @ 3 levels of detail
plot(wavdata, '-');
hold on
plot(wavdata(1:256), 'k-', 'LineWidth', 2);
plot([512.5 512.5],[-2,2],'r-', 'LineWidth', 2);
plot([256.5 256.5],[-2,2],'r-', 'LineWidth', 2);
yline(-1.0, 'b-.', 'LineWidth', 2); yline(1.0, 'b-.', 'LineWidth', 2); grid on; axis tight
title('2nd Decomposition (j = 2) ')

subplot(6,1,4)
wavdata = dwtr(data, 3, hfilt); %WT of data by 'filt' @ 3 levels of detail
plot(wavdata, '-');
hold on
plot(wavdata(1:128), 'k-', 'LineWidth', 2);
plot([512.5 512.5],[-2,2],'r-', 'LineWidth', 2);
plot([256.5 256.5],[-2,2],'r-', 'LineWidth', 2);
plot([128.5 128.5],[-2,2],'r-', 'LineWidth', 2);
yline(-1.0, 'b-.', 'LineWidth', 2); yline(1.0, 'b-.', 'LineWidth', 2); grid on; axis tight
title('3rd Decomposition (j = 3) ')

subplot(6,1,5)
wavdata = dwtr(data, 4, hfilt); %WT of data by 'filt' @ 3 levels of detail
plot(wavdata, '-');
hold on
plot(wavdata(1:64), 'k-', 'LineWidth', 2);
plot([512.5 512.5],[-2,2],'r-', 'LineWidth', 2);
plot([256.5 256.5],[-2,2],'r-', 'LineWidth', 2);
plot([128.5 128.5],[-2,2],'r-', 'LineWidth', 2);
plot([64.5 64.5],[-2,2],'r-', 'LineWidth', 2);
yline(-1.0, 'b-.', 'LineWidth', 2); yline(1.0, 'b-.', 'LineWidth', 2); grid on; axis tight
title('4th Decomposition (j = 4) ')

subplot(6,1,6)
wavdata = dwtr(data, 5, hfilt); %WT of data by 'filt' @ 3 levels of detail
plot(wavdata, '-');
hold on
plot(wavdata(1:32), 'k-', 'LineWidth', 2);
plot([512.5 512.5],[-max(abs(wavdata)),max(abs(wavdata))],'r-', 'LineWidth', 2);
plot([256.5 256.5],[-max(abs(wavdata)),max(abs(wavdata))],'r-', 'LineWidth', 2);
plot([128.5 128.5],[-max(abs(wavdata)),max(abs(wavdata))],'r-', 'LineWidth', 2);
plot([64.5 64.5],[-max(abs(wavdata)),max(abs(wavdata))],'r-', 'LineWidth', 2);
plot([32.5 32.5],[-max(abs(wavdata)),max(abs(wavdata))],'r-', 'LineWidth', 2);
axis([1 N -1.5 1.5])
yline(-1.0, 'b-.', 'LineWidth', 2); yline(1.0, 'b-.', 'LineWidth', 2); grid on; axis tight
title('5th Decomposition (j = 5) ')

saveas(s, './Figures/SR_WTExample.png')
% subplot(3,1,3)
% 
% lambda = 1;
% swt = wavdata .* (abs(wavdata) > lambda );
% dat = idwtr(swt, 5, hfilt); %IWT of data by 'filt' @ 3 levels of detail
% 
% plot(real(dat))
% axis([1 N -1.5 1.5])
% yline(-1.0, 'b.-'); yline(1.0, 'b.-'); grid on
% title('Doppler Recovered')
% %print -depsc 'C:\Brani\Talks\figs\Dem0a.eps'
% 
