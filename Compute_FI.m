function I = Compute_FI(a, b, sigma, theta)
    % Computes the value of the given function
    % Inputs:
    % a, b: Bounds of the interval
    % sigma: Standard deviation of the normal distribution
    % theta: Mean of the normal distribution
    
    % Ensure sigma is positive
    if sigma <= 0
        error('Sigma must be positive');
    end
    
    % Standardize variables for normal distribution
    z_a = (a - theta) / sigma;
    z_b = (b - theta) / sigma;
    
    % PDF and CDF of standard normal distribution
    phi = @(z) (1 / sqrt(2 * pi)) * exp(-0.5 * z.^2);
    Phi = @(z) 0.5 * (1 + erf(z / sqrt(2)));
    
    % Compute the terms of the equation
    term1 = (phi(z_a)^2) / Phi(z_a);
    term2 = ((phi(z_b) - phi(z_a))^2) / (Phi(z_b) - Phi(z_a));
    term3 = (phi(z_b)^2) / (1 - Phi(z_b));
    
    % Combine terms to compute I
    I = (1 / sigma^2) * (term1 + term2 + term3);
end
