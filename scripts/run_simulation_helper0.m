function run_simulation_helper0(params)

params = get_final_params_struct_rk4(params);

param = params(1);

HHWUIS_solver_NCexamples_RK4(param);

isdeletebin = 1;
merge_all_parts_binfile_NC(param, isdeletebin);

end


function params = get_final_params_struct_rk4(params)
% function params = get_final_params_struct(params)
% makes the array struct so that it can be used to run sims in
% parallel 

params = struct('G_K', params.G_K,... 
    'G_Na', params.G_Na,...
    'C', params.C,...
    'g_l_k', params.g_l_k, ... 
    'g_l_na', params.g_l_na, ... % 0.0175
    'g_l_cl', params.g_l_cl, ...
    'dslp', params.dslp,... 
    'gmag', params.gmag, ...
    'gamma', params.gamma,...  % 1e-2: convert from m to cm.  0.0445;
    'beta', params.beta,...    % intracellular and extracellular volum ratio
    'alpha', params.alpha,... 
    'E_Cl', params.E_Cl,...
    'E1', params.E1,...
    'E2', params.E2,...
    'O2Buff', params.O2Buff,... %
    'pmax', params.pmax,...
    'gee', params.gee, 'gii', params.gii, ...    
    'th', params.th, ... %    
     'N1spiketimes_fn', params.N1spiketimes_fn,...
     'N2spiketimes_fn', params.N2spiketimes_fn,...
     'OO1_fn', params.OO1_fn,...
     'OO2_fn', params.OO2_fn,...
     'VE1_fn', params.VE1_fn,...
     'VE2_fn', params.VE2_fn,...
     'KO1_fn', params.KO1_fn,...
     'KO2_fn', params.KO2_fn,...
     'NAI1_fn', params.NAI1_fn,...
     'NAI2_fn', params.NAI2_fn,...
     'S1_fn', params.S1_fn,...
     'S2_fn', params.S2_fn,...
     'KBuff', params.KBuff,...
     'buffersize', params.buffersize, 'durtn', params.durtn, ... %      
     'Nx', params.Nx, 'Nw', params.Nw, ... % 
     'InitialConditions', params.InitialConditions,...          
     'lambda1', params.lambda1,...
     'lambda2', params.lambda2, 'epsilon_o', params.epsilon_o, 'brcount', params.brcount,...
     'outdirname', params.outdirname, ...     
     'con_prob', params.con_prob, 'rand_net_ee', params.rand_net_ee,...
     'rand_net_ie', params.rand_net_ie, 'rand_net_ii', params.rand_net_ii,...
     'rand_net_ei', params.rand_net_ei, ...          
     'netsize', params.netsize, ...     
     'tau1', params.tau1, 'tau2', params.tau2);
end


function merge_all_parts_binfile_NC(param, isdeletebin)

% tot_length = length(params);
outdirname = param.outdirname;

combine_mats_in_binfile_NC(fullfile(outdirname, [param.N1spiketimes_fn, '.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.N2spiketimes_fn, '.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.OO1_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.OO2_fn, '.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.OO1_fn, '.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.OO2_fn, '.std.bin']), isdeletebin);


combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.Isyn1.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.Isyn2.bin']), isdeletebin);


combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.Isyn1.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.Isyn2.std.bin']), isdeletebin);


combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.mge1.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.mge1.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.mgi1.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.mgi1.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.mge2.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.mge2.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.mgi2.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.mgi2.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.KO1_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.KO2_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.NAI1_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.NAI2_fn, '.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.KO1_fn, '.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.KO2_fn, '.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.NAI1_fn, '.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.NAI2_fn, '.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.S1_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.S2_fn, '.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.S1_fn, '.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.S2_fn, '.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.EK1.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.EK2.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.EK1.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.EK2.std.bin']), isdeletebin);

combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.ENa1.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.ENa2.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE1_fn, '.ENa1.std.bin']), isdeletebin);
combine_mats_in_binfile_NC(fullfile(outdirname, [param.VE2_fn, '.ENa2.std.bin']), isdeletebin);

convert_endvals_bin2mat_NC(fullfile(outdirname, [param.VE1_fn, '.endVals.bin']), isdeletebin);

end


function combine_mats_in_binfile_NC(binfile, isdeletebin)
% read all matrices from a .bin file
% ans saves it as .mat file in -v 7.3 format

[fid, msg] = fopen(binfile, 'r');


if fid<1, error(['Error: ', msg, 'File: ' binfile ]), end 


mf_ob = matfile([binfile(1:end-4), '.mat'], 'Writable', true);
mf_ob.mat = [];
sz = fread(fid, 2, 'double').';

fxd_dim2 = sz(2); %second  dim should same always

var_dim1_ed = 0; % variable first dim endpoint

while (~isempty(sz))
    
    if (sz(2) ~= fxd_dim2)
        error('error: dim2 should always be same');
    end
    
    if(sz(1))
        var_dim1_st = var_dim1_ed + 1;
        var_dim1_ed = var_dim1_ed + sz(1);

        mat = fread(fid, sz, 'double');

        mf_ob.mat(var_dim1_st:var_dim1_ed, 1:fxd_dim2) = mat;
    end
    
    sz = fread(fid, 2, 'double').';            
end

fclose(fid);


if isdeletebin
    delete(binfile);
end

end


function convert_endvals_bin2mat_NC(binfile, isdeletebin)
% read all matrices from a .bin file
% ans saves it as .mat file in -v 7.3 format

fid = fopen(binfile, 'r');

mf_ob = matfile([binfile(1:end-4), '.mat'], 'Writable', true);

sz = fread(fid, 2, 'double').';
mf_ob.VE1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.ME1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.NE1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.HE1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.KO1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.NAI1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.OO1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.S1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.X1end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.VE2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.ME2end = fread(fid, sz, 'double');


sz = fread(fid, 2, 'double').';
mf_ob.NE2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.HE2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.KO2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.NAI2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.OO2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.S2end = fread(fid, sz, 'double');

sz = fread(fid, 2, 'double').';
mf_ob.X2end = fread(fid, sz, 'double');

fclose(fid);

if isdeletebin
    delete(binfile);
end


end

