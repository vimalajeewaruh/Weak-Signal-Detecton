clc; clear all; close all force
addpath '/Users/hvimalajeewa2/Documents/Tamu_Documents/TAMU/Resonance/Matlab_codes/NewCodes/'

rng(1,'twister');

N = 1000;
% Lower threshold
sigma = 1.000; 
% Upper threshold
theta = 2.00; 

% Constant Signal
a = linspace(-1, 1.99, N);
% Noise Level
b = linspace(2.01,5, N);

% Collect FI 
Rs = zeros((N));

for j = 1:N    
    for i = 1:N 
       FI = Compute_FI(a(i), b(j), sigma, theta);
        Rs(j,i) = FI;
    end
end 

[l,m] = find(Rs == max(Rs(:)));

[max(Rs(:)) b(l) a(m)]

%% 
figureHandle = figure;

set(figureHandle, 'Units', 'inches');
set(figureHandle, 'Position', [1, 1, 24, 8]);  % 6 inches wide and 4 inches tall

subplot(121)
[xx,yy] = meshgrid(a,b);
X =xx; Y = yy; Z = Rs;

mesh(X, Y, Z);         % Plot the 3D surface
hold on;

% Initialize vectors for storing maximum points
numCols = size(Z, 2);  % Number of columns in Z
xMax = zeros(1, numCols); % To store x-coordinates of maxima
yMax = zeros(1, numCols); % To store y-coordinates of maxima
zMax = zeros(1, numCols); % To store maximum z-values

% Find maximum Z values for each column
for col = 1:numCols
    [zMax(col), rowIdx] = max(Z(:, col)); % Find max Z in the column
    xMax(col) = X(rowIdx, col);           % Corresponding X coordinate
    yMax(col) = Y(rowIdx, col);           % Corresponding Y coordinate
end

% Plot the maximum points
plot3(xMax, yMax, zMax, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Connect the points with a line
plot3(xMax, yMax, zMax, 'b-', 'LineWidth', 2);

hold off;

%xlabel('Signal(\theta)'); ylabel('Amount of Noise (\sigma)'); 
xlabel('Lower threshold(a)'); ylabel('Upper threshold (b)'); 
zlabel('Fisher Information')
title(' Constant signal (\theta = 2) and  Noise level \sigma = 1')

%print -depsc '/Users/dixon/Documents/TAMU/Resonance/Matlab_codes/Figures/Double_thresh_loss.eps'
%print -dpng '/Users/hvimalajeewa2/Documents/Tamu_Documents/TAMU/Resonance/Matlab_codes/Figures/Double_thresh_loss.png'
%
%colormap(pink)    % change color map
%shading interp  

% Lower threshold
a = -2.00; 
% Upper threshold
b = 2.00; 

% Constant Signal
theta = linspace(-2.0, 2.0, N);
% Noise Level
sigma = linspace(.5, 1.5, N);

% Collect FI 
Rs = zeros((N));

for j = 1:N    
    for i = 1:N 
       FI = Compute_FI(a, b, sigma(i), theta(j));
        Rs(j,i) = FI;
    end
end 

[l,m] = find(Rs == max(Rs(:)));

[max(Rs(:)) theta(m) sigma(l)]

%% 

subplot(122)
[xx,yy] = meshgrid(sigma,theta);

X =xx; Y = yy; Z = Rs;
mesh(X, Y, Z);         % Plot the 3D surface
hold on;

% Initialize vectors for storing maximum points
numCols = size(Z, 2);  % Number of columns in Z
xMax = zeros(1, numCols); % To store x-coordinates of maxima
yMax = zeros(1, numCols); % To store y-coordinates of maxima
zMax = zeros(1, numCols); % To store maximum z-values

% Find maximum Z values for each column
for col = 1:numCols
    [zMax(col), rowIdx] = max(Z( col,:)); % Find max Z in the column
    xMax(col) = X(col, rowIdx);           % Corresponding X coordinate
    yMax(col) = Y(col, rowIdx);           % Corresponding Y coordinate
end

% Plot the maximum points
plot3(xMax, yMax, zMax, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

% Connect the points with a line
plot3(xMax, yMax, zMax, 'b-', 'LineWidth', 2);

hold off;

ylabel('Signal(\theta)'); xlabel('Amount of Noise (\sigma)'); 
%xlabel('Lower threshold(a)'); ylabel('Upper threshold (b)'); 
zlabel('Fisher Information')
title(' Lower threshold(a) = -2 and Upper threshold (b) = 2')

%saveas(figureHandle, '/Users/hvimalajeewa2/Documents/Tamu_Documents/TAMU/Resonance/Matlab_codes/NewCodes/Figures/Test_Double_thesh_FI.png')