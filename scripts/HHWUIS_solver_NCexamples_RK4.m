function [VE1, ME1, HE1, NE1, KO1, NAI1, OO1, S1, X1, ISYN1, VE2, ME2, HE2, NE2, KO2, NAI2, OO2, S2, X2, ISYN2, t, options] = HHWUIS_solver_NCexamples_RK4(options) 


    ve1 = squeeze(options.InitialConditions.VE1(end, :)).';
    me1 = squeeze(options.InitialConditions.ME1(end, :)).';
    he1 = squeeze(options.InitialConditions.HE1(end, :)).';
    ne1 = squeeze(options.InitialConditions.NE1(end, :)).';
    ko1 = squeeze(options.InitialConditions.KO1(end, :)).';
    nai1 = squeeze(options.InitialConditions.NAI1(end, :)).';
    oo1 = squeeze(options.InitialConditions.OO1(end, :)).';
    s1 = squeeze(options.InitialConditions.S1(end, :)).';
    x1 = squeeze(options.InitialConditions.X1(end, :)).';
    
    ve2 = squeeze(options.InitialConditions.VE2(end, :)).';
    me2 = squeeze(options.InitialConditions.ME2(end, :)).';
    he2 = squeeze(options.InitialConditions.HE2(end, :)).';
    ne2 = squeeze(options.InitialConditions.NE2(end, :)).';
    ko2 = squeeze(options.InitialConditions.KO2(end, :)).';
    nai2 = squeeze(options.InitialConditions.NAI2(end, :)).';
    oo2 = squeeze(options.InitialConditions.OO2(end, :)).';
    s2 = squeeze(options.InitialConditions.S2(end, :)).';
    x2 = squeeze(options.InitialConditions.X2(end, :)).';


sdt = 0.05;
dt = 10;
samplerate = 1/sdt*1000;
iters = samplerate*options.durtn;

buffsz_samples = options.buffersize*samplerate;


VE1   = zeros(buffsz_samples, options.Nx); %
ME1   = zeros(buffsz_samples, options.Nx); %
HE1   = zeros(buffsz_samples, options.Nx); %
NE1  = zeros(buffsz_samples, options.Nx); %
KO1  = zeros(buffsz_samples, options.Nx); %
NAI1 = zeros(buffsz_samples, options.Nx); 
OO1 = zeros(buffsz_samples, options.Nx); 
S1 = zeros(buffsz_samples, options.Nx); 
X1 = zeros(buffsz_samples, options.Nx); 
ISYN1 = zeros(buffsz_samples, options.Nx); 

VE2   = zeros(buffsz_samples, options.Nw); %
ME2   = zeros(buffsz_samples, options.Nw); %
HE2   = zeros(buffsz_samples, options.Nw); %
NE2  = zeros(buffsz_samples, options.Nw); %
KO2  = zeros(buffsz_samples, options.Nw); %
NAI2 = zeros(buffsz_samples, options.Nw); 
OO2 = zeros(buffsz_samples, options.Nw); 
S2 = zeros(buffsz_samples, options.Nw); 
X2 = zeros(buffsz_samples, options.Nw); 
ISYN2 = zeros(buffsz_samples, options.Nw); 



MGE1 = zeros(buffsz_samples, options.Nx);
MGI1 = zeros(buffsz_samples, options.Nx);
MGE2 = zeros(buffsz_samples, options.Nw);
MGI2 = zeros(buffsz_samples, options.Nw);

N1spike_times = zeros(buffsz_samples/dt, 1);
cellnumber1 = zeros(buffsz_samples/dt, 1);
N2spike_times = zeros(buffsz_samples/dt, 1);
cellnumber2 = zeros(buffsz_samples/dt, 1);
st_and_nn1 = zeros(buffsz_samples/dt, 2);
st_and_nn2 = zeros(buffsz_samples/dt, 2);

Isyn1allN = zeros(buffsz_samples/dt, options.Nx);
Isyn2allN = zeros(buffsz_samples/dt, options.Nw);

VE1_allN = zeros(buffsz_samples/dt, options.Nx);
VE2_allN = zeros(buffsz_samples/dt, options.Nw);

VE1end = zeros(1, options.Nx);
ME1end = zeros(1, options.Nx);
NE1end = zeros(1, options.Nx);
HE1end = zeros(1, options.Nx);
KO1end = zeros(1, options.Nx);
NAI1end = zeros(1, options.Nx);
OO1end = zeros(1, options.Nx);
S1end = zeros(1, options.Nx);
X1end = zeros(1, options.Nx);

VE2end = zeros(1, options.Nw);
ME2end = zeros(1, options.Nw);
NE2end = zeros(1, options.Nw);
HE2end = zeros(1, options.Nw);
KO2end = zeros(1, options.Nw);
NAI2end = zeros(1, options.Nw);
OO2end = zeros(1, options.Nw);
S2end = zeros(1, options.Nw);
X2end = zeros(1, options.Nw);


pmax = options.pmax;
dslp = options.dslp;
KBuff = options.KBuff;
gmag =options.gmag;
beta = options.beta;
G_Na = options.G_Na;
g_l_na = options.g_l_na;
G_K = options.G_K;
g_l_k = options.g_l_k;
con_prob = options.con_prob;
gee = options.gee;
gii = options.gii;
rand_net_ee = options.rand_net_ee;
rand_net_ie = options.rand_net_ie;
rand_net_ii = options.rand_net_ii;
rand_net_ei = options.rand_net_ei;
Nx = options.Nx;
Nw = options.Nw;
E1 = options.E1;
E2 = options.E2;
g_l_cl = options.g_l_cl;
E_Cl = options.E_Cl;
C = options.C;
gamma = options.gamma;
alpha = options.alpha;
lambda1 = options.lambda1;
epsilon_o = options.epsilon_o;
O2Buff = options.O2Buff;
lambda2 = options.lambda2;
tau1 = options.tau1;
tau2 = options.tau2;



sdtt = sdt/2;
sdtt6 = sdt/6;
bcount = 1;
brcount = options.brcount; % buffer reach count to tell how many times buffer has been reached.


for k = 1:iters-0
    if mod(k,10000) == 0
%         display(['Iteration step ' sprintf('%f',k) ' of ' sprintf('%f', iters) ]);
        fprintf(1, 'Iteration step %0.0f of %0.0f\n', k, iters);
    end
    
    
    % first call
    [Fve10, Fme10, Fhe10, Fne10, Fko10, Fnai10, Foo10, Fs10, Fx10, FIsyn10, Fve20, Fme20, Fhe20, Fne20, Fko20, Fnai20, Foo20,Fs20, Fx20, FIsyn20,mge1, mgi1,mge2, mgi2] = HHWUIS_NC_dxdt(ve1,me1,he1,ne1,ko1,nai1,oo1,s1,x1,ve2,me2,he2,ne2,ko2,nai2,oo2,s2,x2, pmax, dslp, KBuff, gmag, beta, G_Na, g_l_na, G_K, g_l_k, con_prob, gee, gii, rand_net_ee, rand_net_ie, rand_net_ii, rand_net_ei, Nx, Nw, E1, E2, g_l_cl, E_Cl, C, gamma, alpha, lambda1, epsilon_o, O2Buff, lambda2, tau1, tau2);         
    
    ve11 = ve1 + (Fve10) * sdtt; 
    me11 = me1 + (Fme10) * sdtt; 
    he11 = he1 + (Fhe10) * sdtt; 
    ne11 = ne1 + (Fne10) * sdtt; 
    ko11 = ko1 + (Fko10) * sdtt; 
    nai11 = nai1 + (Fnai10) * sdtt;
    oo11 = oo1 + (Foo10) * sdtt;
    s11 = s1 + (Fs10) * sdtt;
    x11 = x1 + (Fx10) * sdtt;
        
    ve21 = ve2 + (Fve20) * sdtt; 
    me21 = me2 + (Fme20) * sdtt; 
    he21 = he2 + (Fhe20) * sdtt; 
    ne21 = ne2 + (Fne20) * sdtt; 
    ko21 = ko2 + (Fko20) * sdtt; 
    nai21 = nai2 + (Fnai20) * sdtt; 
    oo21 = oo2 + (Foo20) * sdtt; 
    s21 = s2 + (Fs20) * sdtt; 
    x21 = x2 + (Fx20) * sdtt;     
    
    
    % second call
    [Fve11, Fme11, Fhe11, Fne11, Fko11, Fnai11, Foo11, Fs11, Fx11,~,  Fve21, Fme21, Fhe21, Fne21, Fko21, Fnai21, Foo21, Fs21, Fx21, ~, ~, ~, ~, ~] = HHWUIS_NC_dxdt(ve11,me11,he11,ne11,ko11,nai11,oo11,s11,x11,ve21,me21,he21,ne21,ko21,nai21,oo21,s21,x21, pmax, dslp, KBuff, gmag, beta, G_Na, g_l_na, G_K, g_l_k, con_prob, gee, gii, rand_net_ee, rand_net_ie, rand_net_ii, rand_net_ei, Nx, Nw, E1, E2, g_l_cl, E_Cl, C, gamma, alpha, lambda1, epsilon_o, O2Buff, lambda2, tau1, tau2);
    
    ve12 = ve1 + (Fve11) * sdtt; 
    me12 = me1 + (Fme11) * sdtt; 
    he12 = he1 + (Fhe11) * sdtt; 
    ne12 = ne1 + (Fne11) * sdtt; 
    ko12 = ko1 + (Fko11) * sdtt; 
    nai12 = nai1 + (Fnai11) * sdtt;
    oo12 = oo1 + (Foo11) * sdtt;
    s12 = s1 + (Fs11) * sdtt;
    x12 = x1 + (Fx11) * sdtt;
        
    ve22 = ve2 + (Fve21) * sdtt; 
    me22 = me2 + (Fme21) * sdtt; 
    he22 = he2 + (Fhe21) * sdtt; 
    ne22 = ne2 + (Fne21) * sdtt; 
    ko22 = ko2 + (Fko21) * sdtt; 
    nai22 = nai2 + (Fnai21) * sdtt; 
    oo22 = oo2 + (Foo21) * sdtt; 
    s22 = s2 + (Fs21) * sdtt; 
    x22 = x2 + (Fx21) * sdtt; 
           
    
    
    % third call
    [Fve12, Fme12, Fhe12, Fne12, Fko12, Fnai12, Foo12, Fs12, Fx12,~,  Fve22, Fme22, Fhe22, Fne22, Fko22, Fnai22, Foo22, Fs22, Fx22, ~, ~, ~, ~, ~] = HHWUIS_NC_dxdt(ve12,me12,he12,ne12,ko12,nai12,oo12,s12,x12,ve22,me22,he22,ne22,ko22,nai22,oo22,s22,x22, pmax, dslp, KBuff, gmag, beta, G_Na, g_l_na, G_K, g_l_k, con_prob, gee, gii, rand_net_ee, rand_net_ie, rand_net_ii, rand_net_ei, Nx, Nw, E1, E2, g_l_cl, E_Cl, C, gamma, alpha, lambda1, epsilon_o, O2Buff, lambda2, tau1, tau2);
    
    ve13 = ve1 + (Fve12) * sdt; 
    me13 = me1 + (Fme12) * sdt; 
    he13 = he1 + (Fhe12) * sdt; 
    ne13 = ne1 + (Fne12) * sdt; 
    ko13 = ko1 + (Fko12) * sdt; 
    nai13 = nai1 + (Fnai12) * sdt;
    oo13 = oo1 + (Foo12) * sdt;
    s13 = s1 + (Fs12) * sdt;
    x13 = x1 + (Fx12) * sdt;
        
    ve23 = ve2 + (Fve22) * sdt; 
    me23 = me2 + (Fme22) * sdt; 
    he23 = he2 + (Fhe22) * sdt; 
    ne23 = ne2 + (Fne22) * sdt; 
    ko23 = ko2 + (Fko22) * sdt; 
    nai23 = nai2 + (Fnai22) * sdt; 
    oo23 = oo2 + (Foo22) * sdt; 
    s23 = s2 + (Fs22) * sdt; 
    x23 = x2 + (Fx22) * sdt; 
    
    
    % Forth call
    [Fve13, Fme13, Fhe13, Fne13, Fko13, Fnai13, Foo13, Fs13, Fx13,~,  Fve23, Fme23, Fhe23, Fne23, Fko23, Fnai23, Foo23, Fs23, Fx23, ~, ~, ~, ~, ~] = HHWUIS_NC_dxdt(ve13,me13,he13,ne13,ko13,nai13,oo13,s13,x13,ve23,me23,he23,ne23,ko23,nai23,oo23,s23,x23, pmax, dslp, KBuff, gmag, beta, G_Na, g_l_na, G_K, g_l_k, con_prob, gee, gii, rand_net_ee, rand_net_ie, rand_net_ii, rand_net_ei, Nx, Nw, E1, E2, g_l_cl, E_Cl, C, gamma, alpha, lambda1, epsilon_o, O2Buff, lambda2, tau1, tau2);
    
    nve1 = ve1 + sdtt6 * (Fve10 + 2*Fve11 + 2*Fve12 + Fve13); 
    nme1 = me1 + sdtt6 * (Fme10 + 2*Fme11 + 2*Fme12 + Fme13); 
    nhe1 = he1 + sdtt6 * (Fhe10 + 2*Fhe11 + 2*Fhe12 + Fhe13); 
    nne1 = ne1 + sdtt6 * (Fne10 + 2*Fne11 + 2*Fne12 + Fne13); 
    nko1 = ko1 + sdtt6 * (Fko10 + 2*Fko11 + 2*Fko12 + Fko13); 
    nnai1 = nai1 + sdtt6 * (Fnai10 + 2*Fnai11 + 2*Fnai12 + Fnai13); 
    noo1 = oo1 + sdtt6 * (Foo10 + 2*Foo11 + 2*Foo12 + Foo13); 
    ns1 = s1 + sdtt6 * (Fs10 + 2*Fs11 + 2*Fs12 + Fs13);
    nx1 = x1 + sdtt6 * (Fx10 + 2*Fx11 + 2*Fx12 + Fx13); 
    
    nve2 = ve2 + sdtt6 * (Fve20 + 2*Fve21 + 2*Fve22 + Fve23); 
    nme2 = me2 + sdtt6 * (Fme20 + 2*Fme21 + 2*Fme22 + Fme23); 
    nhe2 = he2 + sdtt6 * (Fhe20 + 2*Fhe21 + 2*Fhe22 + Fhe23); 
    nne2 = ne2 + sdtt6 * (Fne20 + 2*Fne21 + 2*Fne22 + Fne23); 
    nko2 = ko2 + sdtt6 * (Fko20 + 2*Fko21 + 2*Fko22 + Fko23); 
    nnai2 = nai2 + sdtt6 * (Fnai20 + 2*Fnai21 + 2*Fnai22 + Fnai23); 
    noo2 = oo2 + sdtt6 * (Foo20 + 2*Foo21 + 2*Foo22 + Foo23); 
    ns2 = s2 + sdtt6 * (Fs20 + 2*Fs21 + 2*Fs22 + Fs23);
    nx2 = x2 + sdtt6 * (Fx20 + 2*Fx21 + 2*Fx22 + Fx23);             
    

    % removing last dim for coder
    VE1(bcount, :) = nve1;
    ME1(bcount, :) = nme1;
    NE1(bcount, :) = nne1;
    HE1(bcount, :) = nhe1;
    KO1(bcount, :) = nko1;
    NAI1(bcount, :) = nnai1;
    OO1(bcount, :) = noo1;
    S1(bcount, :) = ns1;
    X1(bcount, :) = nx1;
    ISYN1(bcount, :) = FIsyn10;
    MGE1(bcount, :) = mge1;
    MGI1(bcount, :) = mgi1;
    
    
    
    
    VE2(bcount, :) = nve2;
    ME2(bcount, :) = nme2;
    NE2(bcount, :) = nne2;
    HE2(bcount, :) = nhe2;
    KO2(bcount, :) = nko2;
    NAI2(bcount, :) = nnai2;
    OO2(bcount, :) = noo2;
    S2(bcount, :) = ns2;
    X2(bcount, :) = nx2;
    ISYN2(bcount, :) = FIsyn20;
    MGE2(bcount, :) = mge2;
    MGI2(bcount, :) = mgi2;
    

    
    if bcount == buffsz_samples       

       % saving spiketimes -- 
       % TODO: try without downsampling. If there is not a significant decrease in time. keep the non-downsampled one. 
       [N1spike_times, cellnumber1] = get_spiketimes_cellnumber_NC(VE1(1:dt:buffsz_samples, :), options.th, brcount*(buffsz_samples/dt));
             
       if ~isempty(cellnumber1)
        st_and_nn1 = [N1spike_times, cellnumber1];                        
        fid = fopen([options.outdirname,'/',[options.N1spiketimes_fn,'.bin']], 'a');
        fwrite(fid, size(st_and_nn1), 'double');
        fwrite(fid, st_and_nn1, 'double');
        fclose(fid);
       else
        fid = fopen([options.outdirname,'/',[options.N1spiketimes_fn,'.bin']], 'a');
        fwrite(fid, [0,2], 'double');
        fclose(fid);        
       end
       
       [N2spike_times, cellnumber2] = get_spiketimes_cellnumber_NC(VE2(1:dt:buffsz_samples, :), options.th, brcount*(buffsz_samples/dt));
       
       if ~isempty(cellnumber2)
           st_and_nn2 = [N2spike_times, cellnumber2];
           fid = fopen([options.outdirname,'/',[options.N2spiketimes_fn,'.bin']], 'a');
           fwrite(fid, size(st_and_nn2), 'double');
           fwrite(fid, st_and_nn2, 'double');
           fclose(fid);
       else
           fid = fopen([options.outdirname,'/',[options.N2spiketimes_fn,'.bin']], 'a');
           fwrite(fid, [0,2], 'double');
           fclose(fid);
       end

       OO1_ds = mean(OO1(1:dt:buffsz_samples, :),2);
       OO2_ds = mean(OO2(1:dt:buffsz_samples, :),2);
       
       OO1_ds_std = std(OO1(1:dt:buffsz_samples, :),0,2);
       OO2_ds_std = std(OO2(1:dt:buffsz_samples, :),0,2);
       
                                  
       fid = fopen([options.outdirname,'/',[options.OO1_fn,'.bin']], 'a');
       fwrite(fid, size(OO1_ds), 'double');
       fwrite(fid, OO1_ds, 'double');
       fclose(fid);       

       fid = fopen([options.outdirname,'/',[options.OO2_fn,'.bin']], 'a');
       fwrite(fid, size(OO2_ds), 'double');
       fwrite(fid, OO2_ds, 'double');
       fclose(fid);
                     
       
       fid = fopen([options.outdirname,'/',[options.OO1_fn,'.std.bin']], 'a');
       fwrite(fid, size(OO1_ds_std), 'double');
       fwrite(fid, OO1_ds_std, 'double');
       fclose(fid);       

       fid = fopen([options.outdirname,'/',[options.OO2_fn,'.std.bin']], 'a');
       fwrite(fid, size(OO2_ds_std), 'double');
       fwrite(fid, OO2_ds_std, 'double');
       fclose(fid);              
       
       
       % saving EK and ENa
       my_E_K1 = 26.64 * log(KO1./(140.0 + (18.0 - NAI1)));
       my_E_Na1 = 26.64 * log((144.0 - beta*(NAI1 - 18.0))./NAI1);
       
       my_E_K2 = 26.64 * log(KO2./(140.0 + (18.0 - NAI2)));
       my_E_Na2 = 26.64 * log((144.0 - beta*(NAI2 - 18.0))./NAI2);
       
       my_E_K1_ds_mean = mean(my_E_K1(1:dt:buffsz_samples, :),2);
       my_E_Na1_ds_mean = mean(my_E_Na1(1:dt:buffsz_samples, :),2);       
       my_E_K1_ds_std = std(my_E_K1(1:dt:buffsz_samples, :),0,2);
       my_E_Na1_ds_std = std(my_E_Na1(1:dt:buffsz_samples, :),0,2);       
            
       my_E_K2_ds_mean = mean(my_E_K2(1:dt:buffsz_samples, :),2);
       my_E_Na2_ds_mean = mean(my_E_Na2(1:dt:buffsz_samples, :),2);       
       my_E_K2_ds_std = std(my_E_K2(1:dt:buffsz_samples, :),0,2);
       my_E_Na2_ds_std = std(my_E_Na2(1:dt:buffsz_samples, :),0,2);              
       
       fid = fopen([options.outdirname,'/',[options.VE1_fn,'.EK1.bin']], 'a');
       fwrite(fid, size(my_E_K1_ds_mean), 'double');
       fwrite(fid, my_E_K1_ds_mean, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.VE1_fn,'.EK1.std.bin']], 'a');
       fwrite(fid, size(my_E_K1_ds_std), 'double');
       fwrite(fid, my_E_K1_ds_std, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.VE2_fn,'.EK2.bin']], 'a');
       fwrite(fid, size(my_E_K2_ds_mean), 'double');
       fwrite(fid, my_E_K2_ds_mean, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.VE2_fn,'.EK2.std.bin']], 'a');
       fwrite(fid, size(my_E_K2_ds_std), 'double');
       fwrite(fid, my_E_K2_ds_std, 'double');
       fclose(fid);
       
       
       
       fid = fopen([options.outdirname,'/',[options.VE1_fn,'.ENa1.bin']], 'a');
       fwrite(fid, size(my_E_Na1_ds_mean), 'double');
       fwrite(fid, my_E_Na1_ds_mean, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.VE1_fn,'.ENa1.std.bin']], 'a');
       fwrite(fid, size(my_E_Na1_ds_std), 'double');
       fwrite(fid, my_E_Na1_ds_std, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.VE2_fn,'.ENa2.bin']], 'a');
       fwrite(fid, size(my_E_Na2_ds_mean), 'double');
       fwrite(fid, my_E_Na2_ds_mean, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.VE2_fn,'.ENa2.std.bin']], 'a');
       fwrite(fid, size(my_E_Na2_ds_std), 'double');
       fwrite(fid, my_E_Na2_ds_std, 'double');
       fclose(fid);
       
       % other savings
       VE1_ds = mean(VE1(1:dt:buffsz_samples, :),2);
       VE2_ds = mean(VE2(1:dt:buffsz_samples, :),2);
       
       
       VE1_ds_std = std(VE1(1:dt:buffsz_samples, :), 0, 2);
       VE2_ds_std = std(VE2(1:dt:buffsz_samples, :), 0, 2);
       
       
       fid = fopen([options.outdirname,'/',[options.VE1_fn,'.bin']], 'a');
       fwrite(fid, size(VE1_ds), 'double');
       fwrite(fid, VE1_ds, 'double');
       fclose(fid);              

       fid = fopen([options.outdirname,'/',[options.VE2_fn,'.bin']], 'a');
       fwrite(fid, size(VE2_ds), 'double');
       fwrite(fid, VE2_ds, 'double');
       fclose(fid);       
       
       
       fid = fopen([options.outdirname,'/',[options.VE1_fn,'.std.bin']], 'a');
       fwrite(fid, size(VE1_ds_std), 'double');
       fwrite(fid, VE1_ds_std, 'double');
       fclose(fid);              

       fid = fopen([options.outdirname,'/',[options.VE2_fn,'.std.bin']], 'a');
       fwrite(fid, size(VE2_ds_std), 'double');
       fwrite(fid, VE2_ds_std, 'double');
       fclose(fid);       
       
              
        Isyn1 = mean(-ISYN1(1:dt:buffsz_samples, :),2);
        
        Isyn1_std = std(-ISYN1(1:dt:buffsz_samples, :),0, 2);
                
        fid = fopen([options.outdirname,'/',[options.VE1_fn,'.Isyn1.bin']], 'a');
        fwrite(fid, size(Isyn1), 'double');
        fwrite(fid, Isyn1, 'double');
        fclose(fid);
        
        fid = fopen([options.outdirname,'/',[options.VE1_fn,'.Isyn1.std.bin']], 'a');
        fwrite(fid, size(Isyn1_std), 'double');
        fwrite(fid, Isyn1_std, 'double');
        fclose(fid);
        
        Isyn2 = mean(-ISYN2(1:dt:buffsz_samples, :),2);                        
        Isyn2_std = std(-ISYN2(1:dt:buffsz_samples, :),0, 2);
                        
        fid = fopen([options.outdirname,'/',[options.VE2_fn,'.Isyn2.bin']], 'a');
        fwrite(fid, size(Isyn2), 'double');
        fwrite(fid, Isyn2, 'double');
        fclose(fid);             
        
        fid = fopen([options.outdirname,'/',[options.VE2_fn,'.Isyn2.std.bin']], 'a');
        fwrite(fid, size(Isyn2_std), 'double');
        fwrite(fid, Isyn2_std, 'double');
        fclose(fid);         
                        
        mge1 = mean(MGE1(1:dt:buffsz_samples, :),2);
        mge1_std = std(MGE1(1:dt:buffsz_samples, :),0,2);
        
        fid = fopen([options.outdirname,'/',[options.VE1_fn,'.mge1.bin']], 'a');
        fwrite(fid, size(mge1), 'double');
        fwrite(fid, mge1, 'double');
        fclose(fid);
                    
        
        fid = fopen([options.outdirname,'/',[options.VE1_fn,'.mge1.std.bin']], 'a');
        fwrite(fid, size(mge1_std), 'double');
        fwrite(fid, mge1_std, 'double');
        fclose(fid);         
        
        
        mgi1 = mean(MGI1(1:dt:buffsz_samples, :),2);
        mgi1_std = std(MGI1(1:dt:buffsz_samples, :),0,2);
        
        fid = fopen([options.outdirname,'/',[options.VE1_fn,'.mgi1.bin']], 'a');
        fwrite(fid, size(mgi1), 'double');
        fwrite(fid, mgi1, 'double');
        fclose(fid);
                    
        
        fid = fopen([options.outdirname,'/',[options.VE1_fn,'.mgi1.std.bin']], 'a');
        fwrite(fid, size(mgi1_std), 'double');
        fwrite(fid, mgi1_std, 'double');
        fclose(fid);
        
      
        
        mge2 = mean(MGE2(1:dt:buffsz_samples, :),2);
        mge2_std = std(MGE2(1:dt:buffsz_samples, :),0,2);
        
        fid = fopen([options.outdirname,'/',[options.VE2_fn,'.mge2.bin']], 'a');
        fwrite(fid, size(mge2), 'double');
        fwrite(fid, mge2, 'double');
        fclose(fid);
                    
        
        fid = fopen([options.outdirname,'/',[options.VE2_fn,'.mge2.std.bin']], 'a');
        fwrite(fid, size(mge2_std), 'double');
        fwrite(fid, mge2_std, 'double');
        fclose(fid);         
        
        
        mgi2 = mean(MGI2(1:dt:buffsz_samples, :),2);
        mgi2_std = std(MGI2(1:dt:buffsz_samples, :),0,2);
        
        fid = fopen([options.outdirname,'/',[options.VE2_fn,'.mgi2.bin']], 'a');
        fwrite(fid, size(mgi2), 'double');
        fwrite(fid, mgi2, 'double');
        fclose(fid);
                    
        
        fid = fopen([options.outdirname,'/',[options.VE2_fn,'.mgi2.std.bin']], 'a');
        fwrite(fid, size(mgi2_std), 'double');
        fwrite(fid, mgi2_std, 'double');
        fclose(fid);
      
                              
       KO1_ds = mean(KO1(1:dt:buffsz_samples, :),2);
       KO2_ds = mean(KO2(1:dt:buffsz_samples, :),2);
       
       KO1_ds_std = std(KO1(1:dt:buffsz_samples, :), 0, 2);
       KO2_ds_std = std(KO2(1:dt:buffsz_samples, :), 0, 2);
       
             

       fid = fopen([options.outdirname,'/',[options.KO1_fn,'.bin']], 'a');
       fwrite(fid, size(KO1_ds), 'double');
       fwrite(fid, KO1_ds, 'double');
       fclose(fid);       

       fid = fopen([options.outdirname,'/',[options.KO2_fn,'.bin']], 'a');
       fwrite(fid, size(KO2_ds), 'double');
       fwrite(fid, KO2_ds, 'double');
       fclose(fid);     
       
       
       fid = fopen([options.outdirname,'/',[options.KO1_fn,'.std.bin']], 'a');
       fwrite(fid, size(KO1_ds_std), 'double');
       fwrite(fid, KO1_ds_std, 'double');
       fclose(fid);       

       fid = fopen([options.outdirname,'/',[options.KO2_fn,'.std.bin']], 'a');
       fwrite(fid, size(KO2_ds_std), 'double');
       fwrite(fid, KO2_ds_std, 'double');
       fclose(fid);
       
       NAI1_ds = mean(NAI1(1:dt:buffsz_samples, :),2);
       NAI2_ds = mean(NAI2(1:dt:buffsz_samples, :),2);
       
       NAI1_ds_std = std(NAI1(1:dt:buffsz_samples, :), 0, 2);
       NAI2_ds_std = std(NAI2(1:dt:buffsz_samples, :), 0, 2);
       
       fid = fopen([options.outdirname,'/',[options.NAI1_fn,'.bin']], 'a');
       fwrite(fid, size(NAI1_ds), 'double');
       fwrite(fid, NAI1_ds, 'double');
       fclose(fid);

       fid = fopen([options.outdirname,'/',[options.NAI2_fn,'.bin']], 'a');
       fwrite(fid, size(NAI2_ds), 'double');
       fwrite(fid, NAI2_ds, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.NAI1_fn,'.std.bin']], 'a');
       fwrite(fid, size(NAI1_ds_std), 'double');
       fwrite(fid, NAI1_ds_std, 'double');
       fclose(fid);

       fid = fopen([options.outdirname,'/',[options.NAI2_fn,'.std.bin']], 'a');
       fwrite(fid, size(NAI2_ds_std), 'double');
       fwrite(fid, NAI2_ds_std, 'double');
       fclose(fid);
       
       S1_ds = mean(S1(1:dt:buffsz_samples, :),2);
       S2_ds = mean(S2(1:dt:buffsz_samples, :),2);
       
       S1_ds_std = std(S1(1:dt:buffsz_samples, :), 0, 2);
       S2_ds_std = std(S2(1:dt:buffsz_samples, :), 0, 2);                    
       

       fid = fopen([options.outdirname,'/',[options.S1_fn,'.bin']], 'a');
       fwrite(fid, size(S1_ds), 'double');
       fwrite(fid, S1_ds, 'double');
       fclose(fid);

       fid = fopen([options.outdirname,'/',[options.S2_fn,'.bin']], 'a');
       fwrite(fid, size(S2_ds), 'double');
       fwrite(fid, S2_ds, 'double');
       fclose(fid);
       
       fid = fopen([options.outdirname,'/',[options.S1_fn,'.std.bin']], 'a');
       fwrite(fid, size(S1_ds_std), 'double');
       fwrite(fid, S1_ds_std, 'double');
       fclose(fid);

       fid = fopen([options.outdirname,'/',[options.S2_fn,'.std.bin']], 'a');
       fwrite(fid, size(S2_ds_std), 'double');
       fwrite(fid, S2_ds_std, 'double');
       fclose(fid);
       
%        
       bcount = 1;
       brcount = brcount + 1; % increasing the buffer reach count
    else
       bcount = bcount+1;
    end

    

    ve1 = nve1;
    me1 = nme1;
    ne1 = nne1;
    he1 = nhe1;
    ko1 = nko1;
    nai1 = nnai1;
    oo1 = noo1;
    s1 = ns1; %
    x1 = nx1;
    
    ve2 = nve2;
    me2 = nme2;
    ne2 = nne2;
    he2 = nhe2;
    ko2 = nko2;
    nai2 = nnai2;
    oo2 = noo2;
    s2 = ns2; %
    x2 = nx2;
    
end

    VE1end = VE1(end, :);
    ME1end = ME1(end, :);
    NE1end = NE1(end, :);
    HE1end = HE1(end, :);
    KO1end = KO1(end, :);
    NAI1end = NAI1(end, :);
    OO1end = OO1(end, :);
    S1end = S1(end, :);
    X1end = X1(end, :);

    VE2end = VE2(end, :);
    ME2end = ME2(end, :);
    NE2end = NE2(end, :);
    HE2end = HE2(end, :);
    KO2end = KO2(end, :);
    NAI2end = NAI2(end, :);
    OO2end = OO2(end, :);
    S2end = S2(end, :);
    X2end = X2(end, :);
    

fid = fopen([options.outdirname,'/',[options.VE1_fn,'.endVals.bin']], 'w');
fwrite(fid, size(VE1end), 'double');
fwrite(fid, VE1end, 'double');

fwrite(fid, size(ME1end), 'double');
fwrite(fid, ME1end, 'double');

fwrite(fid, size(NE1end), 'double');
fwrite(fid, NE1end, 'double');

fwrite(fid, size(HE1end), 'double');
fwrite(fid, HE1end, 'double');

fwrite(fid, size(KO1end), 'double');
fwrite(fid, KO1end, 'double');

fwrite(fid, size(NAI1end), 'double');
fwrite(fid, NAI1end, 'double');

fwrite(fid, size(OO1end), 'double');
fwrite(fid, OO1end, 'double');

fwrite(fid, size(S1end), 'double');
fwrite(fid, S1end, 'double');

fwrite(fid, size(X1end), 'double');
fwrite(fid, X1end, 'double');



fwrite(fid, size(VE2end), 'double');
fwrite(fid, VE2end, 'double');

fwrite(fid, size(ME2end), 'double');
fwrite(fid, ME2end, 'double');

fwrite(fid, size(NE2end), 'double');
fwrite(fid, NE2end, 'double');

fwrite(fid, size(HE2end), 'double');
fwrite(fid, HE2end, 'double');

fwrite(fid, size(KO2end), 'double');
fwrite(fid, KO2end, 'double');

fwrite(fid, size(NAI2end), 'double');
fwrite(fid, NAI2end, 'double');

fwrite(fid, size(OO2end), 'double');
fwrite(fid, OO2end, 'double');

fwrite(fid, size(S2end), 'double');
fwrite(fid, S2end, 'double');

fwrite(fid, size(X2end), 'double');
fwrite(fid, X2end, 'double');

fclose(fid);

end
                                             
function [Fve1, Fme1, Fhe1, Fne1, Fko1, Fnai1, Foo1, Fs1, Fx1, FIsyn1, Fve2, Fme2, Fhe2, Fne2, Fko2, Fnai2, Foo2, Fs2, Fx2, FIsyn2, mge1, mgi1,mge2, mgi2 ] = HHWUIS_NC_dxdt(ve1,me1,he1,ne1,ko1,nai1,oo1,s1,x1,ve2,me2,he2,ne2,ko2,nai2,oo2,s2,x2, pmax, dslp, KBuff, gmag, beta, G_Na, g_l_na, G_K, g_l_k, con_prob, gee, gii, rand_net_ee, rand_net_ie, rand_net_ii, rand_net_ei, Nx, Nw, E1, E2, g_l_cl, E_Cl, C, gamma, alpha, lambda1, epsilon_o, O2Buff,  lambda2, tau1, tau2)    %#codegen
    
    % pump, glia and diffusion (mM/s)
    p1 = (pmax./(1+exp(-(oo1-20)/3)));
    I_pump1 = (p1./(1.0+exp((25.0-nai1)/3.0))).*(1.0./(1.0+exp(5.5-ko1)));
    Igliapump1 = (p1./(1.0+exp((25.0-18)/3.0)))/3.*(1.0./(1.0+exp(5.5-ko1)));
    I_diff1 = (dslp)*(ko1-KBuff) + ((gmag)./(1.0 + exp((18.0-ko1)/2.5)))+2*Igliapump1;
    
    p2 = (pmax./(1+exp(-(oo2-20)/3)));
    I_pump2 = (p2./(1.0+exp((25.0-nai2)/3.0))).*(1.0./(1.0+exp(5.5-ko2)));
    Igliapump2 = (p2./(1.0+exp((25.0-18)/3.0)))/3 .* (1.0./(1.0+exp(5.5-ko2)));
    I_diff2 = (dslp)*(ko2-KBuff) + ((gmag)./(1.0 + exp((18.0-ko2)/2.5)))+2*Igliapump2;
        
    E_K1 = 26.64 * log(ko1./(140.0 + (18.0 - nai1)));
    E_Na1 = 26.64 * log((144.0 - beta*(nai1 - 18.0))./nai1);
    
    E_K2 = 26.64 * log(ko2./(140.0 + (18.0 - nai2)));
    E_Na2 = 26.64 * log((144.0 - beta*(nai2 - 18.0))./nai2);
        
    % optimized
    INa_E1 = G_Na * (me1.*me1.*me1) .* he1 .* (ve1 - E_Na1) + (g_l_na * (ve1 - E_Na1));
    IK_E1 = G_K * (ne1.*ne1.*ne1.*ne1) .* (ve1 - E_K1) + (g_l_k * (ve1 - E_K1));       

    % optimized
    INa_E2 = G_Na * (me2.*me2.*me2) .* he2 .* (ve2 - E_Na2) + (g_l_na * (ve2 - E_Na2));
    IK_E2 = G_K * (ne2.*ne2.*ne2.*ne2) .* (ve2 - E_K2) + (g_l_k * (ve2 - E_K2));   
    

    x1temp = x1; 
    x2temp = x2;    
    ie_temp = (gee*(s1.*exp(-x1temp/5))); % Nx*1
    ii_temp = (gii*(s2.*exp(-x2temp/5))); % Nw*1

    %%% E to E
    iee = (rand_net_ee * ie_temp);% nx x 1

    %%% I to E
    iie = (rand_net_ie * ii_temp);% nx x 1

    %%% i to i
    iii = (rand_net_ii * ii_temp);% nw x 1

    %%% E to I
    iei = (rand_net_ei * ie_temp);% nw x 1

    Isyn1 = (ve1 - E2) .* iie + (ve1 - E1) .* iee; % nx x 1
    Isyn2 = (ve2 - E1) .* iei + (ve2 - E2) .* iii; % nw x 1
    
    mge1 = ((ve1 - E1) .* iee);
    mgi1 = ((ve1 - E2) .* iie);
    mge2 = ((ve2 - E1) .* iei);
    mgi2 = ((ve2 - E2) .* iii);
            
    
    Fve1 = (- INa_E1 - IK_E1 - (g_l_cl * (ve1 - E_Cl)) - (Isyn1))/C;
    Fme1 = (0.32*(54 + ve1) ./ (1-exp(-(ve1 + 54)/4))) .* (1-me1) - (0.28*(ve1 + 27) ./ (exp((ve1 + 27)/5) - 1)).*me1;
    Fhe1 = (0.128*exp(-(50+ve1)/18)) .* (1-he1) - (4./(1+exp(-(ve1+27)/5))).*he1;
    Fne1 = (0.032*(ve1+52) ./ (1-exp(-(ve1+52)/5))) .* (1-ne1) - (0.5*exp(-(ve1+57)/40)).*ne1;
    Fko1 = 0.001*(gamma*beta*IK_E1 - beta*2.0*I_pump1 - I_diff1);
    Fnai1 = 0.001*(-1.0*gamma*INa_E1 - 3.0*I_pump1);
    Foo1 = 0.001*(-alpha*lambda1*(I_pump1 + Igliapump1) + epsilon_o*(O2Buff-oo1));
    
    %tau1 ; 
    nu1 = zeros(size(ve1));
    nu1((ve1 >= -30) & (ve1 <= -10)) = 0.4;
    Fs1 = (1/tau1)*((20./(1 + exp(-(ve1 + 20)/3))).*(1-s1)-s1);
    Fx1 = nu1.*(ve1 + 50) - 0.4*x1;
    
    Fve2 = (- INa_E2 - IK_E2 - (g_l_cl * (ve2 - E_Cl)) - (Isyn2))/C;
    Fme2 = (0.32*(54 + ve2) ./ (1-exp(-(ve2 + 54)/4))) .* (1-me2) - ( 0.28*(ve2 + 27) ./ (exp((ve2 + 27)/5) - 1)).*me2;
    Fhe2 = (0.128*exp(-(50+ve2)/18)) .* (1-he2) - (4./(1+exp(-(ve2+27)/5))).*he2;
    Fne2 = (0.032*(ve2+52) ./ (1-exp(-(ve2+52)/5))) .* (1-ne2) - (0.5*exp(-(ve2+57)/40)).*ne2;
    Fko2 = 0.001*(gamma*beta*IK_E2 - beta*2.0*I_pump2 - I_diff2);
    Fnai2 = 0.001*(-1.0*gamma*INa_E2 - 3.0*I_pump2);
    Foo2 = 0.001*(-alpha*lambda2*(I_pump2 + Igliapump2) + epsilon_o*(O2Buff-oo2));
    
    %tau2 ;
    nu2 = zeros(size(ve2));
    nu2((ve2 >= -30) & (ve2 <= -10)) = 0.4;
    Fs2 = (1/tau2)*((20./(1 + exp(-(ve2 + 20)/3))).*(1-s2)-s2);
    Fx2 = nu2.*(ve2 + 50) - 0.4*x2;
    
    FIsyn1 = Isyn1;
    FIsyn2 = Isyn2;
    
end


function [N1spike_times, cellnumber] = get_spiketimes_cellnumber_NC(data, th, start_time) %#codegen        
    
    frth = th;
    data_df = data>frth;
    data_df = [zeros(1,size(data_df,2));diff(data_df,1,1)];
    data_spk = zeros(size(data_df));
    data_spk(data_df==1) = 1;

    [N1spike_times, cellnumber] = find(data_spk);
    N1spike_times = N1spike_times + start_time; 

end
