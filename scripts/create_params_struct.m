function params = create_params_struct(regime_type, durtn, outdirpath)
            
    switch regime_type
        case 'AI'            
            params.outdirname = fullfile(outdirpath,'AI');
            params.KBuff{1} = 3.5;
            params.O2Buff{1} = 32;
            ICfilename = './ICs/AI_IC.mat';
        case 'Iso'
            params.outdirname = fullfile(outdirpath,'Iso');
            params.KBuff{1} = 9.5;
            params.O2Buff{1} = 7;
            ICfilename = './ICs/Iso_IC.mat';
        case 'BS'
            params.outdirname = fullfile(outdirpath,'BS');
            params.KBuff{1} = 20;
            params.O2Buff{1} = 7.049;
            ICfilename = './ICs/BS_IC.mat';
        case 'SZ'
            params.outdirname = fullfile(outdirpath,'SZ');
            params.KBuff{1} = 8;
            params.O2Buff{1} = 11.3333;
            ICfilename = './ICs/SZ_IC.mat';
        otherwise
            error('Error: wrong regime_type! Please manipulate the code for additional functionality or contact the authors');
    end
        
    if ~isfolder(outdirpath)
        mkdir(outdirpath);
    end
    
    % creting the output directory 
    % issues warning if the directory already exists
    mkdir(params.outdirname);
    
    % intialization of some default params
    params = initialize_params(params);                         
        
  
    % total duration of simulation
    params.durtn = durtn; %1; % this will be updated for traj sims
       
    % buffersize to for partial savings of output files
    params.buffersize = 1; % in s
    params.brcount = 0;  % 0 for NC sims
                         
    params.InitialConditions = SetInitialConditions(ICfilename);
                          
    save(fullfile(params.outdirname,'params.mat'), 'params');
end



function InitialConditions = SetInitialConditions(myfilename)      
    xx = load(myfilename);
    
    InitialConditions.VE1 = xx.VE1end;
    InitialConditions.VE2 = xx.VE2end;

    InitialConditions.ME1 = xx.ME1end;
    InitialConditions.HE1 = xx.HE1end;
    InitialConditions.NE1 = xx.NE1end;
    InitialConditions.KO1 = xx.KO1end;
    InitialConditions.NAI1 = xx.NAI1end;
    InitialConditions.OO1 = xx.OO1end;
    InitialConditions.S1 = xx.S1end;
    InitialConditions.X1 = xx.X1end;


    InitialConditions.ME2 = xx.ME2end;
    InitialConditions.HE2 = xx.HE2end;
    InitialConditions.NE2 = xx.NE2end;
    InitialConditions.KO2 = xx.KO2end;
    InitialConditions.NAI2 = xx.NAI2end;
    InitialConditions.OO2 = xx.OO2end;
    InitialConditions.S2 = xx.S2end;
    InitialConditions.X2 = xx.X2end;       
end
