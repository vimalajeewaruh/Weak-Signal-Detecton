function smoothed = kernel_smooth_2d(data, kernel_width)
    % Create Gaussian kernel
    [kx, ky] = meshgrid(-kernel_width:kernel_width, -kernel_width:kernel_width);
    gaussian_kernel = exp(-(kx.^2 + ky.^2) / (2 * kernel_width^2));
    gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:)); % Normalize

    % Apply convolution
    smoothed = conv2(data, gaussian_kernel, 'same');
end
