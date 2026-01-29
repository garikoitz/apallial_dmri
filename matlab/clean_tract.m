function clean_tract(tck_file, output_clean, output_sf, maxDist, maxLen, numNodes, maxIter)
% CLEAN_TRACT Clean a TCK tract file and compute SuperFiber
%
% Usage:
%   clean_tract(tck_file, output_clean, output_sf, maxDist, maxLen, numNodes, maxIter)
%
% Inputs:
%   tck_file     - Path to input .tck file
%   output_clean - Path for cleaned tract output
%   output_sf    - Path for SuperFiber output
%   maxDist      - Max Gaussian distance from tract core (std dev), default=3
%   maxLen       - Max fiber length in std dev from mean, default=3
%   numNodes     - Resample each fiber to this many nodes, default=100
%   maxIter      - Maximum cleaning iterations, default=10
%
% Example:
%   clean_tract('/path/to/tracks.tck', '/path/to/tracks_clean.tck', '/path/to/tracks_SF.tck', 3, 3, 100, 10)

% Add paths to required toolboxes
addpath(genpath('/Users/glerma/soft/gitlab/rtp2-pipeline'));

% Set defaults if not provided
if nargin < 4 || isempty(maxDist), maxDist = 3; end
if nargin < 5 || isempty(maxLen), maxLen = 3; end
if nargin < 6 || isempty(numNodes), numNodes = 100; end
if nargin < 7 || isempty(maxIter), maxIter = 10; end

% Convert string inputs to numbers if needed
if ischar(maxDist) || isstring(maxDist), maxDist = str2double(maxDist); end
if ischar(maxLen) || isstring(maxLen), maxLen = str2double(maxLen); end
if ischar(numNodes) || isstring(numNodes), numNodes = str2double(numNodes); end
if ischar(maxIter) || isstring(maxIter), maxIter = str2double(maxIter); end

fprintf('=== Cleaning tract: %s ===\n', tck_file);
fprintf('Parameters: maxDist=%g, maxLen=%g, numNodes=%d, maxIter=%d\n', ...
        maxDist, maxLen, numNodes, maxIter);

% Load the TCK file
fprintf('Loading tract file...\n');
fg = fgRead(tck_file);
fprintf('Original fibers: %d\n', length(fg.fibers));

% Check if we have enough fibers
if length(fg.fibers) < 10
    fprintf('WARNING: Too few fibers (%d), skipping cleaning\n', length(fg.fibers));
    % Save original as "cleaned"
    AFQ_fgWrite(fg, output_clean, 'tck');
    % Create single fiber SF from mean
    if length(fg.fibers) > 0
        meanFiber = mean(cat(3, fg.fibers{:}), 3);
        fg_SF = dtiNewFiberGroup([fg.name '_SuperFiber'], [0 0 155], [], [], {meanFiber});
        AFQ_fgWrite(fg_SF, output_sf, 'tck');
    end
    return;
end

% Clean the fiber group
fprintf('Cleaning fibers...\n');
M = 'mean';  % Use mean for tract core
count = 1;   % Print results
show = 0;    % No visualization

try
    [fg_clean, keep, SuperFiber] = AFQ_removeFiberOutliers(fg, maxDist, maxLen, numNodes, M, count, maxIter, show);
    fprintf('Cleaned fibers: %d (removed %d)\n', length(fg_clean.fibers), length(fg.fibers) - length(fg_clean.fibers));
catch ME
    fprintf('WARNING: Cleaning failed: %s\n', ME.message);
    fprintf('Saving original tract without cleaning\n');
    fg_clean = fg;
    % Compute simple SuperFiber
    [SuperFiber, ~] = dtiComputeSuperFiberRepresentation(fg, [], numNodes);
end

% Create SuperFiber fiber group
fg_SF = dtiNewFiberGroup([fg.name '_SuperFiber'], [0 0 155], [], [], SuperFiber.fibers);

% Save outputs
fprintf('Saving cleaned tract: %s\n', output_clean);
AFQ_fgWrite(fg_clean, output_clean, 'tck');

fprintf('Saving SuperFiber: %s\n', output_sf);
AFQ_fgWrite(fg_SF, output_sf, 'tck');

fprintf('=== Cleaning complete ===\n');

end
