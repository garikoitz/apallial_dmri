%% Standalone script to clean a TCK tract and compute SuperFiber
% Add paths to required toolboxes
addpath(genpath('/Users/glerma/soft/gitlab/rtp2-pipeline'));

%% 1. Load the TCK file
repodir = "/Users/glerma/soft/mbrain/";
datadir = fullfile(repodir,"data");
sub = "129-921-R1L1F";
ses = "260116085433";
dwidir = fullfile(datadir,strcat("derivatives/mrtrix_preproc/sub-", sub, "/ses-", ses, "/dwi/tracts"));
tckfilename = "tracks_Isocortex_Striatum_left.tck";
tck_file = fullfile(dwidir, tckfilename);
fg = fgRead(tck_file);

%% 2. Set cleaning parameters
maxDist  = 2;       % Max Gaussian distance from tract core (std dev)
maxLen   = 2;       % Max fiber length (in std dev from mean length)  
numNodes = 100;     % Resample each fiber to this many nodes
M        = 'mean';  % Use 'mean' or 'median' for tract core
count    = 0;       % Print results to console (1=yes, 0=no)
maxIter  = 10;       % Maximum cleaning iterations
show     = 0;       % Show visualization (1=yes, 0=no)

%% 3. Clean the fiber group and get the SuperFiber
[fg_clean, keep, SuperFiber] = AFQ_removeFiberOutliers(fg, maxDist, maxLen, numNodes, M, count, maxIter, show);

%% 4. Create a fiber group from the SuperFiber (single core fiber)
fg_SF = dtiNewFiberGroup([fg.name '_SuperFiber'], [0 0 155], [], [], SuperFiber.fibers);

%% 5. Save outputs
AFQ_fgWrite(fg_clean, ...
            fullfile(dwidir, strrep(tckfilename,'.tck','_clean.tck')), ...
            'tck');
AFQ_fgWrite(fg_SF, ...
            fullfile(dwidir, ...
            strrep(tckfilename,'.tck','_SF.tck')), 'tck');

%% Optional: Just compute SuperFiber without cleaning
% [SuperFiber_only, fg_resampled] = dtiComputeSuperFiberRepresentation(fg, [], 100);
% fg_SF_only = dtiNewFiberGroup('SuperFiber', [0 0 155], [], [], SuperFiber_only.fibers);
% 
% fprintf('Original fibers: %d\n', length(fg.fibers));
% fprintf('Cleaned fibers: %d\n', length(fg_clean.fibers));
% fprintf('SuperFiber nodes: %d\n', size(SuperFiber.fibers{1}, 2));