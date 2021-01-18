%% Example of the Back-And-Forth method

% Parameters
n = 512;         % The size of the grid
maxIters = 200;  % Maximum iterations for the Back-and-Forth method
TOL = 1e-2;      % Tolerance for the Back-and-Forth method
nt  = 40;        % The number of outer iterations
tau = 0.005;     % Time step for JKO scheme
m   = 3;         % m of internal energy
gamma = 0.01;    % gamma of internal energy
verbose  = 1;    % print out logs 
folder = 'data'; % Set the folder directory to save the data

[x,y] = meshgrid(linspace(0,1,n));

% Define initial density
mu = zeros(n);
idx = (x-0.5).^2 + (y-0.5).^2 < 0.2^2;
mu(idx) = 1;
mu = mu / sum(mu(:)) * n*n;

% Define a quadratic potential
V = 1 * (1 + cos(3*pi*x) .* cos(2*pi*y));

% Define an obstacle
obstacle = zeros(n);

% Run BFM!
result = wgfslow(mu, V, obstacle, m, gamma, maxIters, TOL, nt, tau, folder, verbose);


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

fig = figure;
movieName = 'movie.gif';

for i = 0:nt
        file = fopen(sprintf("%s/rho-%04d.dat", folder, i), 'r');
        rho = fread(file, [n n], 'double');
        imagesc(rho);
        axis square
        set(gca,'XTickLabel',[], 'YTickLabel',[])
        frame = getframe(fig);
        im = frame2im(frame);
        [X,cmap] = rgb2ind(im, 256);

        % Write to the GIF file
        if i == 0
                imwrite(X, cmap, movieName, 'gif', 'DelayTime', 0.1, 'Loopcount', inf);
        else
                imwrite(X, cmap, movieName, 'gif', 'DelayTime', 0.1, 'WriteMode', 'append');
        end
end