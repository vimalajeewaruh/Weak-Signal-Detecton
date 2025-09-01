
function [sig, optimal_parameters_table] = Threshoding2DWTestSig(original_image, sigma_values, kernel_width_values, a, b, num_trials, v, wavelet_level, wavelet_filt, ww)
    % Image size
    image_size = size(original_image);       % Size of the 2D image
    
    % Wavelet filter
    W = Wavmat(wavelet_filt, size(original_image,1), wavelet_level, 0);
    
    % Define necessary functions
    F_inv = @(p) sqrt(2) * erfinv(2 * p - 1);  % Inverse of CDF of Gaussian
    f = @(z) (1 / sqrt(2 * pi)) * exp(-z.^2 / 2);  % PDF of Gaussian
    
    % Initialize variables for optimal parameters search
    best_sigma_sub = 0;
    best_kernel_width_sub = 0;
    min_mse_sub = inf;
    best_sig_sub = zeros(image_size);
    
    best_sigma_sup = 0;
    best_kernel_width_sup = 0;
    min_mse_sup = inf;
    best_sig_sup = zeros(image_size);
    
    best_sigma_joint = 0;
    best_kernel_width_joint = 0;
    min_mse_joint = inf;
    best_sig_joint = zeros(image_size);
    
    % Loop through sigma and kernel width
    for sigma = sigma_values
        for kernel_width = kernel_width_values
            % Initialize accumulators for averaged recoveries
            recovery_sub_avg = zeros(image_size);
            recovery_sup_avg = zeros(image_size);
    
            var_sub_avg = zeros(image_size);
            var_sup_avg = zeros(image_size);
    
            % Repeat noise, thresholding, and recovery for multiple trials
            for trial = 1:num_trials
                % Wavelet Transform of Original Image
                %[wavelet_coeffs, wavelet_dims] = wavedec2(original_image, L, wavelet_name);
                wavelet_coeffs = W * original_image * W';
    
                % Add Gaussian noise in wavelet domain
                noisy_wavelet_coeffs = wavelet_coeffs + sigma * randn(size(wavelet_coeffs));
    
                % Reconstruct noisy image from wavelet domain
                %noisy_image = waverec2(noisy_wavelet_coeffs, wavelet_dims, wavelet_name);
    
                % Thresholding
                sub_threshold = double(noisy_wavelet_coeffs < a);
                sup_threshold = double(noisy_wavelet_coeffs > b);
    
                % Kernel smoothing
                p_a = kernel_smooth_2d(sub_threshold, kernel_width);
                p_b = kernel_smooth_2d(sup_threshold, kernel_width);
    
                % Variance calculations
                %v_a = p_a .* (1 - p_a) ./ (numel(original_image) * f(F_inv(p_a)).^2);
                %v_b = p_b .* (1 - p_b) ./ (numel(original_image) * f(F_inv(p_b)).^2);
    
                xx = max(-v, min(v, F_inv(p_a)));
                v_a = p_a .* (1 - p_a) ./ (numel(original_image) * f(xx).^2);
                xx = max(-v, min(v, F_inv(p_b)));
                v_b = p_b .* (1 - p_b) ./ (numel(original_image) * f(xx).^2);
    
                % Recover signals
                %recovery_sub = a - F_inv(p_a);
                %recovery_sup = b - F_inv(1 - p_b);
    
                recovery_sub = a - max(-v, min(v, F_inv(p_a)));
                recovery_sup = b - min(v, max(-v, F_inv(1 - p_b)));
    
                % Accumulate recoveries
                recovery_sub_avg = recovery_sub_avg + recovery_sub;
                recovery_sup_avg = recovery_sup_avg + recovery_sup;
    
    
                var_sub_avg = var_sub_avg + v_a;
                var_sup_avg = var_sup_avg + v_b;
            end
    
            % Average the recoveries
            recovery_sub_avg = recovery_sub_avg / num_trials;
            recovery_sup_avg = recovery_sup_avg / num_trials;
    
            va = var_sub_avg / num_trials; va(find(~va)) = 1e-3;
            vb = var_sup_avg / num_trials; vb(find(~vb)) = 1e-3;
            
            % Transform back to original domain
            recovery_sub_avg1 = W' * recovery_sub_avg * W; 
            recovery_sub_avg1 = smoothdata(recovery_sub_avg1,'gaussian', ww);
    
            recovery_sup_avg2 = W' * recovery_sup_avg * W; 
            recovery_sup_avg2 = smoothdata(recovery_sup_avg2,'gaussian', ww);
            
            % Compute final joint recovery
            final_joint_recovery = recovery_sub_avg .* (vb ./ (va + vb)) + ...
                             recovery_sup_avg .* (va ./ (va + vb));
    
            final_joint_recovery = W' * final_joint_recovery * W; 
            final_joint_recovery = smoothdata(final_joint_recovery,'gaussian', ww);
           
            % Convert joint recovery to wavelet coefficients
            %joint_recovery_wavelet = wavedec2(joint_recovery, L, wavelet_name);
    
            % Transform back to original domain
            %final_joint_recovery = waverec2(final_joint_recovery, wavelet_dims, wavelet_name);
    
            % Compute MSE
            mse_sub = mean((original_image(:) - recovery_sub_avg1(:)).^2);
            mse_sup = mean((original_image(:) - recovery_sup_avg2(:)).^2);
            mse_joint = mean((original_image(:) - final_joint_recovery(:)).^2);
    
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
                best_sig_joint = final_joint_recovery;
            end
        end
    end
    
    % Compute final MSEs
    %final_mse_sub = mean((original_image(:) - best_sig_sub(:)).^2);
    %final_mse_sup = mean((original_image(:) - best_sig_sup(:)).^2);
    %final_mse_joint = mean((original_image(:) - best_sig_joint(:)).^2);
    
    sig.best_sig_sub = best_sig_sub;
    sig.best_sig_sup = best_sig_sup;
    sig.best_sig_joint = best_sig_joint;
    
    % Display optimal parameters
    % Create a table with the optimal parameters
    methods = {'Sub-thresholding'; 'Sup-thresholding'; 'Double-thresholding'};
    optimal_sigma = [best_sigma_sub; best_sigma_sup; best_sigma_joint];
    optimal_kernel_width = [best_kernel_width_sub; best_kernel_width_sup; best_kernel_width_joint];
    minimum_mse = [min_mse_sub; min_mse_sup; min_mse_joint];
    
    % Combine into a MATLAB table
    optimal_parameters_table = table(methods, optimal_sigma, optimal_kernel_width, minimum_mse, ...
        'VariableNames', {'Method', 'Optimal_Sigma', 'Optimal_Kernel_Width', 'Minimum_MSE'});
    
    % Display the table
    %disp(optimal_parameters_table);
end 