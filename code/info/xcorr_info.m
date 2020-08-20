
% correlation function between info and noninfo neurons
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
tau=20;                                                                % choose between [10,20,50]
K=400;

%%

namea={'V1','V4'};
namep={'tar', 'test'};

task=['compute xcorr between info/noninfo', ' K='  sprintf('%1.0i',K), ' tau ='  sprintf('%1.0i',tau)];
display(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/transfer_mc/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/');
addpath('/home/veronika/synced/transfer_result/weight/tags/');       % for tag

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results
      
loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['tag_info_s+c_', sprintf('%1.0i',K)];
load(loadname3)

%% correlation function between signals of plus and minus neurons

nbses=length(w_all);
ncv=size(w_all{1},1);
lags=-K+1:K-1;

r_info=zeros(nbses,2*K-1);

tic
for sess=1:nbses
   
    weights=w_all{sess};
    
    %% get tags
    tag=tag_info{sess};
    anti_tag=(tag==0);
    tags=cat(2,tag,anti_tag);
    
    %%
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y),spike_diff, spike_same,'UniformOutput', false);                                                
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    %% compute LDR for plus and minus neurons in trials, in every permutation cycle 
    
    r_perm=zeros(nperm,2*K-1);
    for perm=1:nperm
        
        xp=cell(2,1);
        for r=1:2
            
            use_tags=double(find(tags(:,r)));
            
            rp=randperm(sum(J));                                                                       % random order of trials
            spikes_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false);                     % random order to all spike trains    
            
            for cv=1:ncv
                spikes_perm{cv}(:,use_tags,:)=spike_one{cv}(:,use_tags,:);                              % neurons with tag info/not info (1=info, 2=notinfo) have regular spike counts
            end
            
            spikes=spikes_perm;
            
            %%
            
            [x_rec] = reconstruct_1c_fun(weights,spikes,kernel);
            xp{r}=x_rec;                                                                             % collect for plus and minus
                                                                        
        end
        
        %% compute correlation function between plus and minus neurons for each session
    
        x=xp{1}; % minus
        y=xp{2}; % plus
    
        [rxy] = correlation_fun(x,y); % correlation between plus and minus signals in trials and cv, averged across trials and cv
        r_perm(perm,:)=rxy;
        
    end
    
    r_info(sess,:)=nanmean(r_perm);
    
end
toc

%% save

if saveres==1
        
    savename=['xcorr_info_tau',sprintf('%1.0i',tau)];
    savefile='/home/veronika/synced/transfer_result/signal/info/';
    save([savefile,savename],'r_info','lags')
    %clear all
    
end
   


