%% Example of the Back-And-Forth method

% Parameters
n = 512;         % The size of the grid
maxIters = 1000;  % Maximum iterations for the Back-and-Forth method
TOL = 1e-3;      % Tolerance for the Back-and-Forth method
nt  = 40;        % The number of outer iterations
tau = 0.002;     % Time step for JKO scheme
verbose  = 1;    % print out logs 
saveData = 1;    % Save the csv data in figures folder
filename = 'data/mu'; % Set the filename


[x,y] = meshgrid(linspace(0,1,n));

% Define initial density
mu = zeros(n);
idx = (x-0.25).^2 + (y-0.25).^2 < 0.1^2;
mu(idx) = 1;
% mu = mu / sum(mu(:)) * n*n;

% Define a quadratic potential
V = 5 * ((x-0.9).^2 + (y-0.9).^2);

% Define an obstacle
obstacle = zeros(n);
idx = (x-0.5).^2 + (y-0.5).^2 < 0.2^2;
obstacle(idx) = 1;

% Run BFM!
% result = wgfslow(mu, V, obstacle, maxIters, TOL, nt, tau, m, gamma, verbose, saveData, filename);
result = wgfinc(mu, V, obstacle, maxIters, TOL, nt, tau, verbose, saveData, filename);


%% Plot

subplot(1,2,1)
imagesc(mu);
title("initial density")
colormap bone
axis square

hold on
subplot(1,2,2)
imagesc(result);
title("final density")
colormap bone
axis square
hold off

saveas(gcf, "./figures/initial-final.png");

% Save figures in figures folder
figure(1)
for i = 0:nt
    subplot(1,1,1);
    filename_full = sprintf("%s-%04d.dat", filename, i);
    fileId = fopen(filename_full);
    im = fread(fileId,[n n],'double');
    imagesc(im);
    colormap bone
    axis square
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
    filename_full = sprintf("./figures/fig-%d.png", i);
    saveas(gcf, filename_full)
end