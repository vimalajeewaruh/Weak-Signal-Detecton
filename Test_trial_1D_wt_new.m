close all; clear all; clc;
rng(1, 'twister');

% Parameters
x = linspace(pi, 8 * pi, 1024);  % 1D domain
original_signal = sin(x);       % 1D sine wave signal

% sigma_values = linspace(.3, 7.0, 20);   % Range of noise levels
% kernel_width_values = linspace(20, 100, 20); % Range of kernel widths
sigma_values = linspace(.1, 3.0, 50);     % Range of noise levels
kernel_width_values = linspace(.1, 1, 50); % Range of kernel widths

a = -2.0;                      % Lower threshold
b = 2.0;                       % Upper threshold
num_trials = 10;               % Number of trials
v = 4.5;                       % Saturation threshold
wavelet_name = 'db4';          % Wavelet family
wavelet_level = 3;             % Decomposition level

% Define necessary functions
F_inv = @(p) sqrt(2) * erfinv(2 * p - 1);  % Inverse of CDF of Gaussian
f = @(z) (1 / sqrt(2 * pi)) * exp(-z.^2 / 2);  % PDF of Gaussian

% Initialize variables for optimal parameter search
best_sigma_sub = 0;
best_kernel_width_sub = 0;
min_mse_sub = inf;
best_sig_sub = zeros(size(original_signal));

best_sigma_sup = 0;
best_kernel_width_sup = 0;
min_mse_sup = inf;
best_sig_sup = zeros(size(original_signal));

best_sigma_joint = 0;
best_kernel_width_joint = 0;
min_mse_joint = inf;
best_sig_joint = zeros(size(original_signal));

% best_sigma = 0;
% best_kernel_width = 0;
% best_sig = 0;
% min_joint_mse = inf;

% Loop through sigma and kernel width
for sigma = sigma_values
    for kernel_width = kernel_width_values
        % Initialize accumulators for averaged recoveries
        recovery_sub_avg = zeros(size(original_signal));
        recovery_sup_avg = zeros(size(original_signal));

        var_sub_avg = zeros(size(original_signal));
        var_sup_avg = zeros(size(original_signal));

        % Repeat noise, thresholding, and recovery for multiple trials
        for trial = 1:num_trials
            % Wavelet Transform of the Original Signal
            [wavelet_coeffs, wavelet_dims] = wavedec(original_signal, wavelet_level, wavelet_name);

            % Add Gaussian noise in the wavelet domain
            noisy_wavelet_coeffs = wavelet_coeffs + sigma * randn(size(wavelet_coeffs));

            % Reconstruct noisy signal from the wavelet domain
            noisy_signal = waverec(noisy_wavelet_coeffs, wavelet_dims, wavelet_name);

            % Thresholding
            sub_threshold = double(noisy_signal < a);
            sup_threshold = double(noisy_signal > b);

            % Kernel smoothing
            p_a = kernel_smooth_1d(sub_threshold, kernel_width);
            p_b = kernel_smooth_1d(sup_threshold, kernel_width);

            % Variance calculations
            v_a = p_a .* (1 - p_a) ./ (numel(original_signal) * f(F_inv(p_a)).^2);
            v_b = p_b .* (1 - p_b) ./ (numel(original_signal) * f(F_inv(p_b)).^2);

            % Recover signals
            recovery_sub = a - max(-v, min(v, F_inv(p_a)));
            recovery_sup = b - min(v, max(-v, F_inv(1 - p_b)));

            % Convert recoveries back to the wavelet domain
            recovery_sub_wavelet = wavedec(recovery_sub, wavelet_level, wavelet_name);
            recovery_sup_wavelet = wavedec(recovery_sup, wavelet_level, wavelet_name);

            % Reconstruct from wavelet domain to original domain
            recovery_sub_reconstructed = waverec(recovery_sub_wavelet, wavelet_dims, wavelet_name);
            recovery_sup_reconstructed = waverec(recovery_sup_wavelet, wavelet_dims, wavelet_name);

            % Accumulate recoveries
            recovery_sub_avg = recovery_sub_avg + recovery_sub_reconstructed;
            recovery_sup_avg = recovery_sup_avg + recovery_sup_reconstructed;

            var_sub_avg = var_sub_avg + v_a;
            var_sup_avg = var_sup_avg + v_b;
        end

        % Average the recoveries
        recovery_sub_avg = recovery_sub_avg / num_trials;
        recovery_sup_avg = recovery_sup_avg / num_trials;

        va = var_sub_avg / num_trials;
        vb = var_sup_avg / num_trials;

        % Compute joint recovery
        joint_recovery = recovery_sub_avg .* (vb ./ (va + vb)) + ...
                         recovery_sup_avg .* (va ./ (va + vb));

        % Compute MSE for sub, sup, and joint recoveries
        mse_sub = mean((original_signal - recovery_sub_avg).^2);
        mse_sup = mean((original_signal - recovery_sup_avg).^2);
        mse_joint = mean((original_signal - joint_recovery).^2);

        % Update best parameters
        % Update best parameters for sub-thresholding
        if mse_sub < min_mse_sub
            min_mse_sub = mse_sub;
            best_sigma_sub = sigma;
            best_kernel_width_sub = kernel_width;
            best_sig_sub = recovery_sub_avg;
        end

        % Update best parameters for sup-thresholding
        if mse_sup < min_mse_sup
            min_mse_sup = mse_sup;
            best_sigma_sup = sigma;
            best_kernel_width_sup = kernel_width;
            best_sig_sup = recovery_sup_avg;
        end

        % Update best parameters for joint recovery
        if mse_joint < min_mse_joint
            min_mse_joint = mse_joint;
            best_sigma_joint = sigma;
            best_kernel_width_joint = kernel_width;
            best_sig_joint = joint_recovery;
        end
    end
end

%%
% Compute final MSEs
    final_mse_sub = mean((original_signal - best_sig_sub).^2);
    final_mse_sup = mean((original_signal - best_sig_sup).^2);
    final_mse_joint = mean((original_signal - best_sig_joint).^2);
% Display optimal parameters
methods = {'Sub-thresholding'; 'Sup-thresholding'; 'Joint Recovery'};
optimal_sigma = [best_sigma_sub; best_sigma_sup; best_sigma_joint];
optimal_kernel_width = [best_kernel_width_sub; best_kernel_width_sup; best_kernel_width_joint];
minimum_mse = [final_mse_sub; final_mse_sup; final_mse_joint];

% Combine into a MATLAB table
optimal_parameters_table = table(methods, optimal_sigma, optimal_kernel_width, minimum_mse, ...
    'VariableNames', {'Method', 'Optimal_Sigma', 'Optimal_Kernel_Width', 'Minimum_MSE'});

% Display the table
disp(optimal_parameters_table);

% % Final recovery with optimal parameters
% sigma = best_sigma;
% kernel_width = best_kernel_width;
% 
% % Initialize accumulators for final averaged recoveries
% recovery_sub_avg = zeros(size(original_signal));
% recovery_sup_avg = zeros(size(original_signal));
% 
% var_sub_avg = zeros(size(original_signal));
% var_sup_avg = zeros(size(original_signal));
% 
% for trial = 1:num_trials
%     % Wavelet Transform of the Original Signal
%     [wavelet_coeffs, wavelet_dims] = wavedec(original_signal, wavelet_level, wavelet_name);
% 
%     % Add Gaussian noise in the wavelet domain
%     noisy_wavelet_coeffs = wavelet_coeffs + sigma * randn(size(wavelet_coeffs));
% 
%     % Reconstruct noisy signal from the wavelet domain
%     noisy_signal = waverec(noisy_wavelet_coeffs, wavelet_dims, wavelet_name);
% 
%     % Thresholding
%     sub_threshold = double(noisy_signal < a);
%     sup_threshold = double(noisy_signal > b);
% 
%     % Kernel smoothing
%     p_a = kernel_smooth_1d(sub_threshold, kernel_width);
%     p_b = kernel_smooth_1d(sup_threshold, kernel_width);
% 
%     % Variance calculations
%     v_a = p_a .* (1 - p_a) ./ (numel(original_signal) * f(F_inv(p_a)).^2);
%     v_b = p_b .* (1 - p_b) ./ (numel(original_signal) * f(F_inv(p_b)).^2);
% 
%     recovery_sub = a - max(-v, min(v, F_inv(p_a)));
%     recovery_sup = b - min(v, max(-v, F_inv(1 - p_b)));
% 
%     % Convert recoveries back to the wavelet domain
%     recovery_sub_wavelet = wavedec(recovery_sub, wavelet_level, wavelet_name);
%     recovery_sup_wavelet = wavedec(recovery_sup, wavelet_level, wavelet_name);
% 
%     % Reconstruct from wavelet domain to original domain
%     recovery_sub_reconstructed = waverec(recovery_sub_wavelet, wavelet_dims, wavelet_name);
%     recovery_sup_reconstructed = waverec(recovery_sup_wavelet, wavelet_dims, wavelet_name);
% 
%     % Accumulate recoveries
%     recovery_sub_avg = recovery_sub_avg + recovery_sub_reconstructed;
%     recovery_sup_avg = recovery_sup_avg + recovery_sup_reconstructed;
% 
%     var_sub_avg = var_sub_avg + v_a;
%     var_sup_avg = var_sup_avg + v_b;
% end
% 
% recovery_sub_avg = recovery_sub_avg / num_trials;
% recovery_sup_avg = recovery_sup_avg / num_trials;
% 
% va = var_sub_avg / num_trials;
% vb = var_sup_avg / num_trials;
% 
% % Compute final joint recovery
% final_joint_recovery = recovery_sub_avg .* (vb ./ (va + vb)) + ...
%                        recovery_sup_avg .* (va ./ (va + vb));
% 
% % Compute final MSEs
% final_mse_sub = mean((original_signal - recovery_sub_avg).^2);
% final_mse_sup = mean((original_signal - recovery_sup_avg).^2);
% final_mse_joint = mean((original_signal - final_joint_recovery).^2);
% 
% % Display final MSEs
% fprintf('Final MSE (Sub-thresholding): %.6f\n', final_mse_sub);
% fprintf('Final MSE (Sup-thresholding): %.6f\n', final_mse_sup);
% fprintf('Final MSE (Joint Recovery): %.6f\n', final_mse_joint);

% Plot results
figure;

% Subplot 1: Original and noisy signals
subplot(2, 1, 1);
plot(x, original_signal, 'b', 'LineWidth', 1.5, 'DisplayName', 'Original Signal');
hold on;
yline(a, 'k--','LineWidth', 1.5, 'DisplayName', 'Lower threshold'); 
yline(b, 'k--','LineWidth', 1.5, 'DisplayName', 'Upper threshold'); 
plot(x, original_signal + best_sigma_joint * randn(size(x)), 'r.', 'LineWidth', 1.2, 'DisplayName', sprintf('Noisy Signal (\\sigma = %.2f)', best_sigma_joint));
title('Original and Noisy Signals');
xlabel('x');
ylabel('Amplitude');
legend('Location', 'best');
grid on; axis tight; ylim([-5 5]);
hold off;

% Subplot 2: Actual and recovered signals
subplot(2, 1, 2);
plot(x, original_signal, 'b', 'LineWidth', 1.5, 'DisplayName', 'Original Signal');
hold on;
plot(x, best_sig_sub, 'g', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Sub-threshold Recovery (\\sigma = %.2f, h = %.2f, MSE = %.6f)', best_sigma_sub, best_kernel_width_sub, min_mse_sub));
plot(x, best_sig_sup, 'm', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Sup-threshold Recovery (\\sigma = %.2f, h = %.2f, MSE = %.6f)', best_sigma_sup, best_kernel_width_sup, min_mse_sup));
plot(x, best_sig_joint, 'k', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Joint Recovery (\\sigma = %.2f, h = %.2f, MSE = %.6f)', best_sigma_joint, best_kernel_width_joint, min_mse_joint));
title('Recovered Signals from Sub, Sup, and Double Thresholding');
xlabel('x');
ylabel('Amplitude');
legend('Location', 'best');
grid on; axis tight;
hold off;

