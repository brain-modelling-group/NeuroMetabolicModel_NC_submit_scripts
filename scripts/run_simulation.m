function run_simulation(regime_type, durtn, outdirpath)
%

params = create_params_struct(regime_type, durtn, outdirpath);

run_simulation_helper0(params);

save_firing_rates(params);

disp('Job Completed... !!Yay!!');

end
