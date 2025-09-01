function smoothed = kernel_smooth_1d(data, kernel_width)
    % Create Gaussian kernel
    kernel_range = -kernel_width:0.01:kernel_width;
    gaussian_kernel = exp(-(kernel_range.^2) / (2 * kernel_width^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel); % Normalize

    % Apply convolution
    smoothed = conv(data, gaussian_kernel, 'same');
end