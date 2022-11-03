addpath('scripts');

% Examples tuned to save the simulations in the AI, Iso, BS, and SZ regimes

% AI
regime_type = 'AI'; % select either AI/Iso/BS/SZ
durtn = 1;    % provide the duration of the simulation
outdirpath = './output_dir'; % provide the output directory to save the simulation output;

tic; % for getting runtime of the simulation ~ Elapsed time is 15.450183 seconds.
run_simulation(regime_type, durtn, outdirpath); % run the simulation;
toc;


%Iso
regime_type = 'Iso'; % select either AI/Iso/BS/SZ
durtn = 1;    % provide the duration of the simulation
outdirpath = './output_dir'; % provide the output directory to save the simulation output;

tic; % for getting runtime of the simulation ~ Elapsed time is 15.334093 seconds.
run_simulation(regime_type, durtn, outdirpath); % run the simulation;
toc;


% BS
regime_type = 'BS'; % select either AI/Iso/BS/SZ
durtn = 50;    % provide the duration of the simulation
outdirpath = './output_dir'; % provide the output directory to save the simulation output;

tic; % for getting runtime of the simulation ~ Elapsed time is 722.761436 seconds.
run_simulation(regime_type, durtn, outdirpath); % run the simulation;
toc;


% SZ
regime_type = 'SZ'; % select either AI/Iso/BS/SZ
durtn = 60;    % provide the duration of the simulation
outdirpath = './output_dir'; % provide the output directory to save the simulation output;

tic; % for getting runtime of the simulation ~ Elapsed time is 719.905010 seconds 
run_simulation(regime_type, durtn, outdirpath); % run the simulation;
toc;

rmpath("scripts/");
