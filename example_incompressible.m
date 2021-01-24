%% Example of the Back-And-Forth method
% 
% For more explanation see also the documentation here:
% https://wasserstein-gradient-flows.netlify.app/
% 


n = 1024;             % The size of the n x n grid
maxIters = 300;      % Maximum number of BFM iterations
TOL = 5e-2;          % Tolerance for BFM
nt  = 150;           % Number of outer iterations
tau = 0.005;         % Time step in the JKO scheme
folder = 'data';     % Output directory
verbose  = 1;        % Print out logs


[x,y] = meshgrid(linspace(0,1,n));

% Initial density
rhoInitial = zeros(n);
idx = (x-0.25).^2 + (y-0.75).^2 < 0.15^2;
rhoInitial(idx) = 1;

% Potential
V = 3 * ((x-0.7).^2 + (y-0.3).^2);

% Add an obstacle
obstacle = zeros(n);
idx = (x-0.5).^2 + (y-0.5).^2 < 0.1^2;
obstacle(idx) = 1;

% Plots
subplot(1,3,1)
contourf(x, y, rhoInitial)
title('Initial density')
axis('square')

subplot(1,3,2)
contourf(x, y, V);
title('Potential V')
axis('square')

subplot(1,3,3)
contourf(x, y, obstacle)
colormap(gca, 'copper')
title('Obstacle')
axis('square')

%% Run BFM!
rhoFinal = wgfinc(rhoInitial, V, obstacle, maxIters, TOL, nt, tau, folder, verbose);


%% Make movie

fig = figure;
movieName = 'movie.gif';

for i = 0:nt
    file = fopen(sprintf("%s/rho-%04d.dat", folder, i), 'r');
    rho = fread(file, [n n], 'double');
    imagesc(rho)
    axis xy square
    set(gca,'XTickLabel',[], 'YTickLabel',[])

    frame = getframe(fig);
    im = frame2im(frame);
    [X,cmap] = rgb2ind(im, 256);

    % Write to the GIF file
    if i == 0
        imwrite(X, cmap, movieName, 'gif', 'Loopcount',inf, 'DelayTime',0.001);
    else
        imwrite(X, cmap, movieName, 'gif', 'WriteMode','append', 'DelayTime',0.001);
    end
end