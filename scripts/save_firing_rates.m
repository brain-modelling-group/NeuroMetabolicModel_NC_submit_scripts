function save_firing_rates(params)

zzz = 1;
    
xx = load(fullfile(params.outdirname, [params.N1spiketimes_fn{zzz}, '.mat']));
st_and_nn1_full = xx.mat;
xx = load(fullfile(params.outdirname, [params.N2spiketimes_fn{zzz}, '.mat']));
st_and_nn2_full = xx.mat;

st = 0;
ed = params.durtn;

binSz = 2; % 1 ms
[N1_totfr, N1_ft] = totalfr_from_spiketimes_cc(st_and_nn1_full, binSz, 1, 2000, ed,st); % in spikes/ms
[N2_totfr, N2_ft] = totalfr_from_spiketimes_cc(st_and_nn2_full, binSz, 1, 2000, ed,st); % in spikes/ms

N1_totfr = N1_totfr * 1e3; % in spikes/s
N2_totfr = N2_totfr * 1e3; % in spikes/s

ft = N1_ft;

save(fullfile(params.outdirname, [params.VE1_fn{zzz}, '.fr', num2str(binSz),'binSz.mat']), 'N1_totfr', 'N2_totfr', 'ft');

end


function [tot_fr, ft] = totalfr_from_spiketimes_cc(st_and_nn, winSize, step, sf, durtn, start_time)

    if isempty(st_and_nn)
       max_t = durtn*sf;
       st_t = (start_time)*sf + 1;
       sts = (1:step:(max_t-st_t+1-winSize+step)).';
       tot_fr = zeros(length(sts), 1); 
       ft = sts/sf + start_time;
       return;
    end
       
    max_t = durtn*sf;
    st_t = (start_time)*sf + 1;
   
    % extracting st > start_time
    st_and_nn = st_and_nn(st_and_nn(:,1)>=st_t,:);
    
        % extracting st > start_time
    st_and_nn = st_and_nn(st_and_nn(:,1)<=max_t,:);
    
    
    % correcting indexs to start from 1
    st_and_nn(:,1) = st_and_nn(:,1) - st_t +1;

   
    tot_stnn = zeros(max_t, 1);
        
    for i = 1 : size(st_and_nn(:, 1))
%         i
        tot_stnn(st_and_nn(i,1)) =  tot_stnn(st_and_nn(i,1)) + 1;      
    end

    % bin wise smoothing
    sts = (1:step:(max_t-st_t+1-winSize+step)).';
    tot_fr = zeros(length(sts), 1);         
    for sti = 1 : length(sts)        
        tot_fr(sti) =  sum(tot_stnn(sts(sti):(sts(sti) + winSize - 1))); % spikes/ms
    end

    ft = sts/sf + start_time;
end    