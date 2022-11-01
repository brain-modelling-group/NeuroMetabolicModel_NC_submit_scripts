function params = initialize_params(params)
params.G_K = 25; 

params.G_Na = 30; 
G_L = 0.05; 

params.C = 1.0;

params.g_l_k = G_L; 

params.g_l_na = 0.35*G_L; % 0.0175

params.g_l_cl = G_L;

params.dslp =0.33; % K diffusion rate
params.gmag=8; % this is G_glia in the eqs--- glia uptake strength of potassium
params.gamma =0.0445;  % 1e-2: convert from m to cm.  0.0445;
params.beta = 7;    % intracellular and extracellular volum ratio

% Oxygen Dynamics related parameters
params.alpha = 5.3; 

% Fixed Reversal Potentials
chloride_oe = 130.0; chloride_ie = 6.0;
params.E_Cl = 26.64*log(chloride_ie/chloride_oe);

% reversal synaptic potential
params.E1 = 0;
params.E2 = -80;

params.pmax = 1.25; %default value

params.lambda1 = 1;
params.lambda2 = 0.5;

params.epsilon_o = 0.17;

% synaptic time constants
params.tau1 = 4;
params.tau2 = 8;


params.gee = 0.022;
params.gii = 0.374;

% netsize and stufff
params.netsize = 400;

% Getting the number of excitatory neurons and the number of inhibitory
% neurons.
params.Nx =  round(params.netsize* 0.8);
params.Nw =  round(params.netsize* (1-0.8));

% this relevant for output/saving
params.iterable = {{'KBuff', 'O2Buff'}};

% here we get the rand- connection matrix
params.myDegree =  80;

params.con_prob = params.myDegree/params.netsize;

rng(9999); % just some rand seed    
% E to E
params.rand_net_ee = sparse(generate_random_conmat(params.Nx, params.Nx, params.con_prob, 1)); % nx x nx
% I to E
params.rand_net_ie = sparse(generate_random_conmat(params.Nx, params.Nw, params.con_prob, 0)); % nx x nw    
% i to i
params.rand_net_ii = sparse(generate_random_conmat(params.Nw, params.Nw, params.con_prob, 1)); % nw x nw
% E to I
params.rand_net_ei = sparse(generate_random_conmat(params.Nw, params.Nx, params.con_prob, 0)); % nw x nx
rng('default'); %5 setting back the rand seed to default;

% to get all the output filenames
params = get_outfilenames(params);

% threshold to consider a spike 
params.th = -40;

end


function con_mat = generate_random_conmat(d1, d2, con_prob, msk_on)
% function con_mat = generate_random_conmat(d1, d2, con_prob, msk_on)
% function to generate random connection matrix bases on connection
% probability

    % 1mconp
    tempcp = 1-con_prob;
    % From E to E
    con_mat = rand(d1, d2);
    con_mat(con_mat>=tempcp) = 1;
    con_mat(con_mat<tempcp) = 0;
    
    % if mask is on then self connections are off or else they are on
    if msk_on
        msk = (ones(d1,d2) -eye(d1,d2));
        con_mat = con_mat.*msk;
    end
    
end


function params = get_outfilenames(params)

KBuff = params.KBuff{1};
O2Buff = params.O2Buff{1};
            
kkk =1;
params.N1spiketimes_fn{kkk,1} = fullfile( ['N1_spike_times_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
params.N2spiketimes_fn{kkk,1} = fullfile( ['N2_spike_times_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);

params.OO1_fn{kkk,1} = fullfile( ['OO1_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
params.OO2_fn{kkk,1} = fullfile( ['OO2_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);


params.VE1_fn{kkk,1} = fullfile( ['VE1_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
params.VE2_fn{kkk,1} = fullfile( ['VE2_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);

params.KO1_fn{kkk,1} = fullfile( ['KO1_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
params.KO2_fn{kkk,1} = fullfile( ['KO2_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);


params.NAI1_fn{kkk,1} = fullfile( ['NAI1_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
params.NAI2_fn{kkk,1} = fullfile( ['NAI2_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);


params.S1_fn{kkk,1} = fullfile( ['S1_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
params.S2_fn{kkk,1} = fullfile( ['S2_KBuff', num2str(KBuff),'_O2Buff', num2str(O2Buff), '.txt']);
            
end



