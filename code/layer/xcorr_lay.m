% correlation function between plus and minus neurons

clear all 
close all
clc

place=1;
if place==1
    saveres=0;
else
    saveres=1;
end

nperm=1000;

ba=1;
period=2;

tau=20;                                                                % choose between [10,20,50]
K=400;

%% names

namea={'V1','V4'};
namelay={'SG','G','IG'};
namep={'tar','test'};

name2lay=cell(3,1);%name 2lay
idx1=[1,1,2];
idx2=[2,3,3];
for i=1:3
    
    lay1=namelay{idx1(i)};
    lay2=namelay{idx2(i)};
    
    name2lay{i}=[lay1,' & ',lay2];
end

task=['compute xcorr for signals in pairs of layers ', ' K='  sprintf('%1.0i',K)];
display(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                          % support
lambda=1/tau;                                                                         % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/transfer_mc/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/synced/transfer_result/signal/layer/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

% load results
loadname=['weight_bac_',namea{ba},'_',namep{period}, '_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_' sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['ncell_lay_',namea{ba}];
load(loadname3);

%% correlation function between signals of plus and minus neurons

nbses=length(w_all);
ncv=size(w_all{1},1);

tic

r_layer=zeros(nbses,3,2*K-1);
lags=-K+1:K-1;

for sess=1:nbses
    
    weights=w_all{sess};
    ncl=ncell_lay(sess,:);
    s=cumsum([0,ncl]);
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y),spike_diff, spike_same,'UniformOutput', false);                                                
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    %% compute signal in layers
    
    r_perm=zeros(nperm,3,2*K-1);
    for perm=1:nperm
        
        xp=cell(3,1);
        for lay = 1:3
            
            rp=randperm(sum(J));                                                                       % random order of trials
            spike_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false); 
            
            delta = s(lay) + 1 : s(lay+1);
            for cv=1:ncv
                spike_perm{cv}(:,delta,:)=spike_one{cv}(:,delta,:);                                   % neurons with sign r (1=neg, 2=pos) have regular spike counts
            end
            
            spikes=spike_perm;
            [x_rec] = reconstruct_1c_fun(weights,spikes,kernel);                               % compute reconstruction in trials
            xp{lay}=x_rec;                                                                   % collect across permutations
            
        end
        
        %% compute correlation function between the population signal in layer 1 and layer 2
        for i=1:3
            x=xp{idx1(i)}; % layer 1
            y=xp{idx2(i)}; % layer 2
            
            [rxy] = correlation_fun(x,y); % correlation between pop signal for a pair oflayers, averged across trials and cv
            r_perm(perm,i,:)=rxy;
            
        end
        
    end
    %%
    r_layer(sess,:,:)=squeeze(nanmean(r_perm));
    
end
toc


    
%% save

if saveres==1
        
    savename=['xcorr_layer_tau',sprintf('%1.0i',tau)];
    savefile='/home/veronika/synced/transfer_result/signal/layer/';
    save([savefile,savename],'r_layer','lags','name2lay')
    %clear all
end
    

%%

