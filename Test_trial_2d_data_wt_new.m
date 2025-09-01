close all; clear all; clc
%rng(1,'twister')
randn('seed',10)
%%
lw = 2.5;  set(0, 'DefaultAxesFontSize', 16);
fs = 15;  msize = 15;

%% Original Images 
% load lena
% original_image = (lena - min(min(lena))  )./(max(max(lena))-min(min(lena)));
image_size = [512, 512];                     % Size of the 2D image
x = linspace(0, 8*pi, image_size(1));
y = linspace(0, 8*pi, image_size(2));
[X, Y] = meshgrid(x, y);
original_image = sin(X) .* cos(Y);           % Example 2D signal

% original_image = (original_image - min(min(original_image))  )./(max(max(original_image))-min(min(original_image)));

%%
% sigma_values = linspace(.25, 2.0, 20);        % Range of noise levels
% kernel_width_values = linspace(.001, 15, 20); % Range of kernel widths  
% num_trials = 100.0;                           % Number of trials
% v  = 6.50;                                    % Saturation threshold
% ww = 70.0;  

sigma_values = linspace(.5,5.50, 20);        % Range of noise levels
kernel_width_values = linspace(.1, 1.50, 20); % Range of kernel widths  
num_trials = 100.0;                           % Number of trials
v  = 6.0;                                    % Saturation threshold
ww = 80.0; 

%% Wavelet family

wtype='Symmlet'; filtersize = 6;
wavelet_filt = MakeONFilter(wtype, filtersize); % making wavelet filters
wavelet_level = 3;

%% Thresholds

% Lower threshold
if min(original_image) <= 0
   a = 02.0*min(original_image(:)); 
else% Lower threshold
   a = -02.0*min(original_image(:));
end

% Upper threshold
b = 02.0*max(original_image(:));  
    
%fprintf('Processing image: %s\n', image_name);
        
[sig_data, optimal_parameters_table_data] = Threshoding2DTestSig(original_image, sigma_values, ...
                                kernel_width_values, a, b, num_trials, v);

% Results for this image
optimal_parameters_table_data

[sig_wt, optimal_parameters_table_wt] = Threshoding2DWTestSig(original_image, sigma_values, ...
                                kernel_width_values, a, b, num_trials, v, wavelet_level, wavelet_filt, ww);

% Store results for this image
optimal_parameters_table_wt
    
%% Plot results    
figureHandle = figure;

set(figureHandle, 'Units', 'inches');
set(figureHandle, 'Position', [1, 1, 24, 8]);  % 6 inches wide and 4 inches tall

subplot(2,5, 1)
colormap summer
imagesc(original_image); % Original signal
title('Original Image', 'FontSize', fs);
        
subplot(2,5, 2)
      
noisy_image = original_image + optimal_parameters_table_data.Optimal_Sigma(3)*randn(size(original_image));
colormap summer
imagesc(noisy_image); % Original signal
title('Noisy Image', 'FontSize', fs);
        
subplot(2,5, 3)
colormap summer
imagesc(sig_data.best_sig_sub); % Original signal
title('Sub-threshold Recovery', 'FontSize', fs);
        
subplot(2,5, 4)
colormap summer
imagesc(sig_data.best_sig_sup); % Original signal
title('Sup-threshold Recovery', 'FontSize', fs);
    
subplot(2,5, 5)
colormap summer
imagesc(sig_data.best_sig_joint); % Original signal
title('Double-threshold Recovery', 'FontSize', fs);

subplot(2,5, 7)

W = Wavmat(wavelet_filt, size(original_image,1), wavelet_level, 0);
wavelet_coeffs = W * original_image * W';
noisy_image = wavelet_coeffs + optimal_parameters_table_wt.Optimal_Sigma(3)*randn(size(original_image));
colormap summer
imagesc(noisy_image); % Original signal
title('Noisy Image', 'FontSize', fs);
        
subplot(2,5, 8)
colormap summer
imagesc(sig_wt.best_sig_sub); % Original signal
title('Sub-threshold Recovery', 'FontSize', fs);
        
subplot(2,5, 9)
colormap summer
imagesc(sig_wt.best_sig_sup); % Original signal
title('Sup-threshold Recovery', 'FontSize', fs);
    
subplot(2,5, 10)
colormap summer
imagesc(sig_wt.best_sig_joint); % Original signal
title('Double-threshold Recovery', 'FontSize', fs);

%saveas(figureHandle, './Figures/Image_data_wt_2d.png')