
% correlation function between plus and minus neurons
% for all trials

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
tau=20;                                                                        % choose between [10,20,50]
T=400;

%%

namea={'V1','V4'};
namep={'tar', 'test'};

task=['compute xcorr between plus and minus neurons ', ' T='  sprintf('%1.0i',T), ' tau ='  sprintf('%1.0i',tau)];
display(task)

% exponential kernel for convolution
L=100;                                                                          % length of the kernel                                  
tau_vec=0:L;                                                                    % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results
      
loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',T)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
load(loadname2,'spikes_inm','spikes_cnm')


%% correlation function between signals of plus and minus neurons

nbses=length(w_all);
ncv=size(w_all{1},1);
lags=-T+1:T-1;

r_sign=zeros(nbses,2*T-1);

tic
parfor sess=1:10%:nbses
   
    weights=w_all{sess};
    
    %% get sign
    tell_sign=cell(2,ncv);
    for cv=1:ncv
        tell_sign{1,cv}=find(weights(cv,:)<0);              % get the index of negative/positive weights 
        tell_sign{2,cv}=find(weights(cv,:)>0);
    end
    %%
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y),spike_diff, spike_same,'UniformOutput', false);                                                
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    %% compute LDR for plus and minus neurons in trials, in every permutation cycle 
    
    r_perm=zeros(nperm,2*T-1);
    for perm=1:nperm
        
        xp=cell(2,1);
        for sgn=1:2
            
            rp=randperm(sum(J));                                                                       % random order of trials
            spikes_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false);                     % random order to all spike trains    
            
            for cv=1:ncv
                idx_sign=tell_sign{sgn,cv};
                spikes_perm{cv}(:,idx_sign,:)=spike_one{cv}(:,idx_sign,:);                              % neurons with sign r (1=neg, 2=pos) have regular spike counts
            end
            
            spikes=spikes_perm;
            
            %%
            
            [x_rec] = reconstruct_1c_fun(weights,spikes,kernel);
            xp{sgn}=x_rec;                                                                             % collect for plus and minus
                                                                        
        end
        
        %% compute correlation function between plus and minus neurons for each session
    
        x=xp{1}; % minus
        y=xp{2}; % plus
    
        [rxy] = correlation_fun(x,y); % correlation between plus and minus signals in trials and cv, averged across trials and cv
        r_perm(perm,:)=rxy;
        
    end
    
    r_sign(sess,:)=nanmean(r_perm);    % mean across permutations
    
end
toc

%% save

if saveres==1
        
    savename=['xcorr_sign_tau',sprintf('%1.0i',tau),'_',sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/signal/correlation/';
    save([savefile,savename],'r_sign','lags')
    clear variables
    
end
    

%%

