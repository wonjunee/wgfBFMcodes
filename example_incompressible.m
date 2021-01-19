%% Example of the Back-And-Forth method

% Parameters
n = 512;         % The size of the grid
maxIters = 300;  % Maximum iterations for the Back-and-Forth method
TOL = 3e-2;      % Tolerance for the Back-and-Forth method
nt  = 200;        % The number of outer iterations
tau = 0.005;     % Time step for JKO scheme
verbose  = 1;    % print out logs 
folder = 'data'; % Set the filename

[x,y] = meshgrid(linspace(0,1,n));

% Define initial density
mu = zeros(n);
idx = (x-0.25).^2 + (y-0.25).^2 < 0.15^2;
mu(idx) = 1;

% Define a quadratic potential
V = 3 * ((x-0.7).^2 + (y-0.7).^2);

% Define an obstacle
obstacle = zeros(n);
idx = (x-0.5).^2 + (y-0.5).^2 < 0.1^2;
obstacle(idx) = 1;


%% Run BFM!

result = wgfinc(mu, V, obstacle, maxIters, TOL, nt, tau, folder, verbose);


%% Plot

subplot(1,2,1)
imagesc(mu);
title("initial density")
% colormap bone
axis square

hold on
subplot(1,2,2)
imagesc(result);
title("final density")
% colormap bone
axis square
hold off

subplot(1,3,1)
contourf(x,y,mu);
title("initial density")
% colormap bone
axis square

hold on
subplot(1,3,2)
contourf(x,y,V);
title("V")
% colormap bone
axis square

hold on
subplot(1,3,3)
contourf(x,y,obstacle);
title("Obstacle")
colormap(gca, 'copper')
axis square

hold off

saveas(gcf, "./figures/initial-final.png");


fig = figure;
movieName = 'movie.gif';

for i = 0:nt
    file = fopen(sprintf("%s/rho-%04d.dat", folder, i), 'r');
    rho = fread(file, [n n], 'double');
    imagesc(rho)
    axis square
    set(gca,'XTickLabel',[], 'YTickLabel',[])

    frame = getframe(fig);
    im = frame2im(frame);
    [X,cmap] = rgb2ind(im, 256);

    % Write to the GIF file
    if i == 0
        imwrite(X, cmap, movieName, 'gif', 'Loopcount',inf, 'DelayTime',0.03);
    else
        imwrite(X, cmap, movieName, 'gif', 'WriteMode','append', 'DelayTime',0.03);
    end
end