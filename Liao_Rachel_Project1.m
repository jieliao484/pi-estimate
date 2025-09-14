
% (1) for loop
rng(123) % set seed for reproducibility
totalpoints=1e5; % set total number of the random points to 1e5
inside=0; % infalut number of the points inside the circle is 0

% Initialize the timer
tic; % begin time count
for i = 1:totalpoints
    x = rand; % generate random x coordinate
    y = rand; % generate random y coordinate
    if x^2 + y^2 <= 1
        inside = inside + 1; % increment count if point is inside the circle
    end
    piEstimate(i) = (inside / i) * 4; % estimate the value of pi
    piErr(i) = abs(piEstimate(i) - pi);
end

elapsedTime = toc; % end time count

fprintf('Single run with N = %d finished in %.3fs. Final pi_hat=%.6f (error=%.3e)\n',...
    totalpoints, elapsedTime, piEstimate(end), piErr(end));

% plots for the estimate pi in single run
subplot(2,3,1);
plot(1:totalpoints, piEstimate, 'LineWidth',0.5); hold on;
yline(pi, '--', 'True \pi', 'LabelVerticalAlignment','bottom')
xlabel('Number of points'); ylabel('\pi estimate');
title('Monte Carlo Estimate of \pi'); grid on;

% Plot the error over the number of points
subplot(2,3,2);
plot(1:totalpoints, piErr, 'LineWidth',0.5);
xlabel('Number of points'); ylabel('Error');
title('Error in \pi Estimate'); grid on;

%%
numbersb = [1e3 1e4 3e5 1e6 1e7];
kmax = numel(numbersb);

piEstimate = zeros(kmax, 1);
piErr = zeros(kmax, 1);
elapsedTime = zeros(kmax, 1);
% Loop over different numbers of points for multiple runs
for j = 1:kmax
    totalpointsb = numbersb(j);
    insideb = 0; % reset count for each totalpoints
    t0 = tic; % begin time count for each run
    for k = 1:totalpointsb
        x = rand; % generate random x coordinate
        y = rand; % generate random y coordinate
        if x^2 + y^2 <= 1
            insideb = insideb + 1; % increment count if point is inside the circle
        end
    end
    piEstimate(j) = (insideb / totalpointsb) * 4; % estimate the value of pi
    piErr(j) = abs(piEstimate(j) - pi);
    elapsedTime(j) = toc(t0); % end time count
    fprintf('N = %-8d  pi_hat = %.6f  |error| = %.3e  time = %.3fs\n', ...
        totalpointsb, piEstimate(j), piErr(j), elapsedTime(j));
end

% Plot the error of estimate pi for different numbers of points
subplot(2,3,4);
plot(numbersb, piErr, 'o-', 'LineWidth', 0.5);
xlabel('Number of points');
ylabel('Error of estimate \pi');
title('Error of estimate \pi vs Number of Points');
grid on;
% Plot the elapsed time for different numbers of points
subplot(2,3,5);
plot(numbersb, elapsedTime, 'o-', 'LineWidth', 0.5);
xlabel('Number of points');
ylabel('Elapsed Time (s)');
title('Elapsed Time vs Number of Points');
grid on;
% Plot the error of estimate pi for different elapsed time
subplot(2,3,6);
plot(elapsedTime, piErr, 'o-', 'LineWidth', 0.5);
xlabel('Elapsed Time (s)');
ylabel('Error of estimate \pi');
title('Elapsed Time vs Error of estimate \pi');
grid on;


%%
%(2) while loop
function piest = user_pi_precision(sigfigs, alpha, maxnumber, batchsize, seed)
% alpha, maxnumber, batchsize and seed are optional to input

% set Defaults if user omitted arguments
if nargin < 2 || isempty(alpha), alpha = 0.05; % 95% CI
end
if nargin < 3 || isempty(maxnumber), maxnumber = 1e7; % maximum loop 1e7 numbers
end
if nargin < 4 || isempty(batchsize), batchsize = 5000; %batch size 5000
end
if nargin >= 5 && ~isempty(seed), rng(seed); % reproducibility
end

% Validate sigfigs
if any(sigfigs < 1 | sigfigs ~= floor(sigfigs))
    error('sigfigs must be positive integers (scalar or vector).');
end
sigfigs = sort(unique(sigfigs)); % sorted vector
S = numel(sigfigs);

% Relative tolerances for significant figures: 0.5*10^(1-s)
relTol = 0.5 .* 10.^(1 - sigfigs);

% z-quantile for two-sided CI
if exist('norminv','file') == 2
    z = norminv(1 - alpha/2);
else
% erfcinv relation: erfc(z/sqrt(2)) = alpha  =>  z = sqrt(2)*erfcinv(alpha)
    z = sqrt(2) * erfcinv(alpha);
end

% Counters, start with zero counts and mark all targets "not yet achieved"
N = 0; % total points 0 so far
inside = 0; % total inside 0 in circle
achievedN = nan(1, S);  % creates a 1-by-S row vector filled with N 
achievedPi = nan(1, S); % the pi estimate at that N
achievedMask = false(1, S);% whether each target is already achieved (true/false)

% Figure setup
figure('Name','Monte Carlo \pi Estimation','Color','w');
clf; hold on; axis([-1 1 -1 1]); axis square;
box on; grid on;
xlabel('x'); ylabel('y');
title('Random Points in a circle');
t = linspace(0, 2*pi, 600);
plot(cos(t), sin(t), 'k-', 'LineWidth', 1.1);

% Prepare animated plotting
Inside = scatter(nan, nan, 8, 'filled', 'MarkerFaceColor', [0 0.45 0.74], ...
                      'MarkerEdgeColor', 'none'); % blue-ish
Outside = scatter(nan, nan, 8, 'filled', 'MarkerFaceColor', [0.85 0.33 0.10], ...
                       'MarkerEdgeColor', 'none'); % orange-ish
legend({'circle boundary','Inside','Outside'}, 'Location','southoutside');

% Text handle for live status (pi estimate + CI half-width)
Text = text(0.02, 0.97, '', 'Units','normalized', 'VerticalAlignment','top', ...
                 'FontName','Consolas', 'FontSize', 10);

% While loop sampling 
tStart = tic;
maxKeep = 200000; % only display the last 200000 points to save memory
X_in = []; Y_in = [];
X_out = []; Y_out = [];

while (~all(achievedMask)) && (N < maxnumber)
    % Generate a batch of random points in [-1,1] x [-1,1]
    b = min(batchsize, maxnumber - N);
    x = 2*rand(b, 1)-1 ; y = 2*rand(b, 1)-1 ;
    in = (x.^2 + y.^2) <= 1;  % inside indicator

    inside = inside + sum(in);
    N = N + b;

    % Update display buffers (keep newest maxKeep points)
    X_in  = [X_in;  x(in)];   Y_in  = [Y_in;  y(in)];
    X_out = [X_out; x(~in)];  Y_out = [Y_out; y(~in)];
    if numel(X_in) > maxKeep
       drop = numel(X_in) - maxKeep;
       X_in(1:drop) = []; Y_in(1:drop) = [];
    end
    if numel(X_out) > maxKeep
       drop = numel(X_out) - maxKeep;
       X_out(1:drop) = []; Y_out(1:drop) = [];
    end

     % Push new data to the scatter plots
     Inside.XData = X_in; Outside.XData = X_out;
     Inside.YData = Y_in; Outside.YData = Y_out;

     % Compute current estimate, standard error, and CI half-width
     p_hat = inside / N; % estimated proportion inside
     pi_hat = 4 * p_hat; % estimate of pi
     se_hat = 4 * sqrt(p_hat * (1 - p_hat) / max(N,1));  % standard error
     ci_half = z * se_hat; % half-width of (1-alpha) CI

     % Check each sigfig target using relative CI half-width
     for sidx = 1:S
         if ~achievedMask(sidx)
            if ci_half <= relTol(sidx) * pi_hat
                 achievedMask(sidx) = true;
                 achievedN(sidx) = N;
                 achievedPi(sidx) = pi_hat;
            end
          end
      end

     % Update the on-plot status text
     msg = sprintf(['N = %d\n\\pi_hat = %.10f\nSE = %.3g   CI_{%.0f%%} half-width = %.3g'], ...
                       N, pi_hat, se_hat, (1-alpha)*100, ci_half);
     set(Text, 'String', msg);
     drawnow limitrate;
end
elapsed = toc(tStart); % end time count
piest = 4 * (inside / N); % Final estimate

fprintf('\n Monte Carlo π Estimation (alpha = %.3f) \n', alpha);
fprintf('Final N = %d   (elapsed %.3f s)\n', N, elapsed);
fprintf('Final pi_hat = %.12f\n', piest);
fprintf('Significant-figure achievements (relative CI half-width criterion):\n');
for sidx = 1:S
    s = sigfigs(sidx);
    tol = relTol(sidx);
    if ~isnan(achievedN(sidx))
      % Format to s significant figures:
        str_sfig = sprintf(['%.' num2str(s) 'g'], achievedPi(sidx));
        fprintf('  %d s.f.: N = %-10d  pi ≈ %s   (tol = %.3g)\n', s, achievedN(sidx), str_sfig, tol);
    else
        fprintf('  %d s.f.: NOT achieved up to maxN = %d (tol = %.3g)\n', s, maxnumber, tol);
    end
end

% Annotate final value (to the highest requested precision) on the plot
s_max = max(sigfigs);
pi_str = sprintf(['%.' num2str(s_max) 'g'], piest);
annot = sprintf('\\pi ≈ %s  (N = %d)', pi_str, N);
text(0.02, 0.05, annot, 'Units','normalized', 'FontSize', 12,'FontWeight','bold');
title(sprintf('Random Points (Reached up to %d significant figure%s)', ...
          s_max, 's'*(s_max~=1) + ''*(s_max==1)));

end

% user only change this part below.
pi_est = user_pi_precision(4);
