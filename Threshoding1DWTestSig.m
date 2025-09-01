
function [sig, optimal_parameters_table] = Threshoding1DWTestSig(original_signal, sigma_values, kernel_width_values, a, b, num_trials, v, wavelet_level, wavelet_name, ww)
    
    % Wavelet matrix
    n = length(original_signal); J = log2(n);  coarsest = J - wavelet_level;
    %W = Wavmat(wavelet_name,n, J-wavelet_level, 0);
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
                %[wavelet_coeffs, wavelet_dims] = wavedec(original_signal, wavelet_level, wavelet_name);
                wavelet_coeffs = dwtr(original_signal, coarsest, wavelet_name);
    
                % Add Gaussian noise in the wavelet domain
                noisy_wavelet_coeffs = wavelet_coeffs + sigma * randn(size(wavelet_coeffs));
    
                % Reconstruct noisy signal from the wavelet domain
                %noisy_signal = waverec(noisy_wavelet_coeffs, wavelet_dims, wavelet_name);
                 noisy_signal = wavelet_coeffs + sigma * randn(size(wavelet_coeffs));
    
                % Thresholding
                sub_threshold = double(noisy_signal < a);
                sup_threshold = double(noisy_signal > b);
    
                % Kernel smoothing
                p_a = kernel_smooth_1d(sub_threshold, kernel_width);
                p_b = kernel_smooth_1d(sup_threshold, kernel_width);
    
                % Variance calculations
                %v_a = p_a .* (1 - p_a) ./ (numel(original_signal) * f(F_inv(p_a)).^2);
                %v_b = p_b .* (1 - p_b) ./ (numel(original_signal) * f(F_inv(p_b)).^2);
                
                xx = max(-v, min(v, F_inv(p_a)));
                 v_a = p_a .* (1 - p_a) ./ (numel(original_signal) * f(xx).^2);
                 xx = max(-v, min(v, F_inv(p_b)));
                 v_b = p_b .* (1 - p_b) ./ (numel(original_signal) * f(xx).^2);

                % Recover signals
                %recovery_sub = a - F_inv(p_a);
                %recovery_sup = b - F_inv(1 - p_b);
    
                recovery_sub = a - max(-v, min(v, F_inv(p_a)));
                recovery_sup = b - min(v, max(-v, F_inv(1 - p_b)));
    
                recovery_sub_avg = recovery_sub_avg + recovery_sub;%_reconstructed;
                recovery_sup_avg = recovery_sup_avg + recovery_sup;%_reconstructed;
    
                var_sub_avg = var_sub_avg + v_a;
                var_sup_avg = var_sup_avg + v_b;
            end
    
            % Average the recoveries
            recovery_sub_avg = recovery_sub_avg / num_trials;
            recovery_sup_avg = recovery_sup_avg / num_trials;
    
            va = var_sub_avg / num_trials;
            vb = var_sup_avg / num_trials;

            %va(find(~va)) = 1e-3; vb(find(~vb)) = 1e-3;
    
            % Compute joint recovery
            %recovery_sub_avg = idwtr(recovery_sub_avg, coarsest, wavelet_name);
            %recovery_sub_avg = smoothdata(recovery_sub_avg,'gaussian', ww);

            %recovery_sup_avg = idwtr(recovery_sup_avg, coarsest, wavelet_name);
            %recovery_sup_avg = smoothdata(recovery_sup_avg,'gaussian', ww);


            %recovery_sub_avg = idwtr(recovery_sub_avg, coarsest, wavelet_name);
            %recovery_sup_avg = idwtr(recovery_sup_avg, coarsest, wavelet_name);
            %joint_recovery   = idwtr(joint_avg, coarsest, wavelet_name);

            joint_recovery = recovery_sub_avg .* (vb ./ (va + vb)) + ...
                             recovery_sup_avg .* (va ./ (va + vb));
            
            joint_recovery = idwtr(joint_recovery, coarsest, wavelet_name);
            joint_recovery = smoothdata(joint_recovery,'gaussian', ww);

            % sub and sup recoveries
            recovery_sub_avg = idwtr(recovery_sub_avg, coarsest, wavelet_name);
            recovery_sub_avg = smoothdata(recovery_sub_avg,'gaussian', ww);

            recovery_sup_avg = idwtr(recovery_sup_avg, coarsest, wavelet_name);
            recovery_sup_avg = smoothdata(recovery_sup_avg,'gaussian', ww);

            
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

    % best_sig_sub = smoothdata(best_sig_sub,'gaussian', ww);
    % best_sig_sup = smoothdata(best_sig_sup,'gaussian', ww);
    % best_sig_joint = smoothdata(best_sig_joint,'gaussian', ww);

    % Compute final MSEs
    final_mse_sub = mean((original_signal - best_sig_sub).^2);
    final_mse_sup = mean((original_signal - best_sig_sup).^2);
    final_mse_joint = mean((original_signal - best_sig_joint).^2);

    sig.best_sig_sub = best_sig_sub;
    sig.best_sig_sup = best_sig_sup;
    sig.best_sig_joint = best_sig_joint;
    
    % Display optimal parameters
    methods = {'Sub-thresholding'; 'Sup-thresholding'; 'Double-thresholding'};
    optimal_sigma = [best_sigma_sub; best_sigma_sup; best_sigma_joint];
    optimal_kernel_width = [best_kernel_width_sub; best_kernel_width_sup; best_kernel_width_joint];
    minimum_mse = [final_mse_sub; final_mse_sup; final_mse_joint];
    
    % Combine into a MATLAB table
    optimal_parameters_table = table(methods, optimal_sigma, optimal_kernel_width, minimum_mse, ...
        'VariableNames', {'Method', 'Optimal_Sigma', 'Optimal_Kernel_Width', 'Minimum_MSE'});
end