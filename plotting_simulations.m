function plotting_simulations
addpath('./scripts');

outdirpath = './output_dir';

%% plotting AI
regime_type = 'AI';
plotting_helper0(regime_type, outdirpath);

%% plotting BS
regime_type = 'BS';
plotting_helper0(regime_type, outdirpath);

%% plotting SZ
regime_type = 'SZ';
plotting_helper0(regime_type, outdirpath);

% %% plotting Iso
% regime_type = 'Iso';
% plotting_helper0(regime_type, outdirpath);

rmpath('./scripts');
end

function plotting_helper0(regime_type, outdirpath)

outdirname = fullfile(outdirpath, regime_type);

iter = 1;
load(fullfile(outdirname, 'params.mat')); 
params.outdirname = outdirname;
st = 0; 
ed = params.durtn; 


figure; 
set(gcf,'units', 'inches');
set(gcf,'paperpositionmode', 'auto');
set(gcf,'Position', [0,0, 8.3, 11.7/3*2]);
set(gcf,'units', 'pixels');
set(gcf,'color','w');

p = panel();

ttls = lower('ABCDEFGHIJKLMNOPQRSTUVWXYZ');

p.pack('v', {2/5, 3/5*1/3, 3/5*1/3, 3/5*1/3});

x = load(fullfile(params.outdirname,[params.N1spiketimes_fn{iter}, '.mat']));
st_and_nn1_full = x.mat;
x = load(fullfile(params.outdirname,[params.N2spiketimes_fn{iter}, '.mat']));
st_and_nn2_full = x.mat;   

ed_t = ed*2000 -1;
st_t = (st)*2000 + 1;

if ~isempty(st_and_nn1_full)
    st_and_nn1 = st_and_nn1_full(st_and_nn1_full(:,1)>=st_t,:);    
    st_and_nn1 = st_and_nn1(st_and_nn1(:,1)<ed_t,:);    
else
    st_and_nn1 = [];
end

if ~isempty(st_and_nn1_full)
    st_and_nn2 = st_and_nn2_full(st_and_nn2_full(:,1)>=st_t,:);
    st_and_nn2 = st_and_nn2(st_and_nn2(:,1)<ed_t,:);
else
    st_and_nn2 = [];
end

[sorted_st_and_nn1, Idx] = sort_rasterspikes(st_and_nn1, params.Nx); 
[sorted_st_and_nn2, Idw] = sort_rasterspikes(st_and_nn2, params.Nw); 



%% Raster
p(1).select();

hold on;
if ~isempty(sorted_st_and_nn1)
msz = 1e-1;
plot(sorted_st_and_nn1(:,1)/2000, sorted_st_and_nn1(:,2), 'o', 'markersize', msz, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(sorted_st_and_nn2(:,1)/2000, sorted_st_and_nn2(:,2)+params.Nx, 'o', 'markersize', msz, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
end

xlim([st, ed]); 
ylim([1,params.Nx+params.Nw]);
ylabel('Neurons');

ttl = ttls(1);
xl = -0.075;
text(xl,1, ttl, 'units', 'normalized', 'fontweight', 'bold', 'verticalalignment', 'bottom');

set(gca,'box','off')

ax = gca;
ax.XAxis.Visible = 'off';


%% F.R.
p(2).select();

binSz = 2; % 1 ms
xx = load(fullfile(params.outdirname,[params.VE1_fn{iter},'.fr', num2str(binSz),'binSz.mat']));


N1_totfr = xx.N1_totfr(st_t:ed_t); %
N2_totfr = xx.N2_totfr(st_t:ed_t);

ft = xx.ft(st_t:ed_t) - xx.ft(st_t);

N1_avgfr = (N1_totfr)/(params.Nx); % in spikes/(s.neuron)
N2_avgfr = (N2_totfr)/(params.Nw); % in spikes/(s.neuron)



plot(ft, N1_avgfr, 'k');
hold on;
plot(ft, -N2_avgfr, 'r');
hold on;
xlim([st, ed]); 

ylabel({'F.R.' ,'{(spikes/s/cell)}'});

xticklabels('');


ttl = ttls(2);
xl = -0.075;
text(xl,1, ttl, 'units', 'normalized', 'fontweight', 'bold', 'verticalalignment', 'bottom');
set(gca,'box','off')


ax = gca;
ax.XAxis.Visible = 'off';



%% P.S.C.
p(3).select();


x = load(fullfile(params.outdirname,[params.VE1_fn{iter}, '.Isyn1.mat']));
Isyn1 = x.mat(st_t:ed_t);


x = load(fullfile(params.outdirname,[params.OO1_fn{iter}, '.mat']));
OO1 = x.mat(st_t:ed_t);
x = load(fullfile(params.outdirname,[params.OO2_fn{iter}, '.mat']));
OO2 = x.mat(st_t:ed_t);



yyaxis left;
plot(((st_t:ed_t))/2000, Isyn1, 'k');
ylabel({'P.S.C.', '{(\muA / cm^2)}'});


yyaxis right;
plot(((st_t:ed_t))/2000, (OO2+OO1)/2, 'linewidth', 2, 'color', [31, 161, 184]/255);
ylabel({'[O_2]_o (mg/L)'});
xlim([st, ed]); 

yyaxis left;

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [31, 161, 184]/255;

ttl = ttls(3);
xl = -0.075;
text(xl,1, ttl, 'units', 'normalized', 'fontweight', 'bold', 'verticalalignment', 'bottom');
set(gca,'box','off')

ax = gca;
ax.XAxis.Visible = 'off';


%% E.I. balance
p(4).select();

x = load(fullfile(params.outdirname,[params.VE1_fn{iter}, '.mge1.mat']));
mge1 = x.mat(st_t:ed_t);
x = load(fullfile(params.outdirname,[params.VE1_fn{iter}, '.mgi1.mat']));
mgi1 = x.mat(st_t:ed_t);

x = load(fullfile(params.outdirname,[params.VE2_fn{iter}, '.mge2.mat']));
mge2 = x.mat(st_t:ed_t);
x = load(fullfile(params.outdirname,[params.VE2_fn{iter}, '.mgi2.mat']));
mgi2 = x.mat(st_t:ed_t);

pp1 = plot((st_t:ed_t)/2000, log10(-mge1./mgi1), 'k', 'linewidth',  3);
hold on
pp2 = plot((st_t:ed_t)/2000, log10(-mge2./mgi2), 'r', 'linewidth',  0.5);
hold off
ylim([-2,2]);

ylabel({'E-I bal.', '{log_{10}(-EPSC / IPSC)}'});
xlabel('Time (s)');


ttl = ttls(4);
xl = -0.075;
text(xl,1, ttl, 'units', 'normalized', 'fontweight', 'bold', 'verticalalignment', 'bottom');
set(gca,'box','off')

ax = gca;
ax.XAxis.Visible = 'on';

p.fontsize = 8;
p(1).marginbottom = 3;
p(2).marginbottom = 3;
p(3).marginbottom = 3;
p(4).marginbottom = 3;
p.margin = [14, 7.7, 16, 5.3];

p.title(regime_type);

end


function [sorted_st_and_nn, I, counts] = sort_rasterspikes(st_and_nn, N)
% st_and_nn: has spikes_times in dim 1 and neuron number in dim 2
% here counts(i) contains the spike counts of ith neuron
if ~isempty(st_and_nn)

    counts = return_counts(st_and_nn(:,2), N);    
    [~, I] = sort(counts);    
    sorted_st_and_nn = sort_standnn(st_and_nn, I);

else
    
    counts = nan;
    I = nan;    
    sorted_st_and_nn = [];
end

end


function counts = return_counts(A, N)
 counts =zeros(N,1);
 for i = 1:N
   counts(i,1) = sum(A==i); % number of times each unique value is repeated
 end
end

function sorted_st_and_nn = sort_standnn(st_and_nn, I)
    sorted_st_and_nn = st_and_nn;
    for i = 1:length(I)
       sorted_st_and_nn((st_and_nn(:,2) == I(i)), 2) = i; 
    end
end
