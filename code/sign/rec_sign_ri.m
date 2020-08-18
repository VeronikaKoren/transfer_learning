
clear all 
close all
%clc

place=1;

if place==1
    saveres=0;
else
    saveres=1;
end

nperm=1000;

ba=1;
period=2;
tau_prime=20;                                                                % choose between [10,20,40]
K=500;

%%

namebeh={'different','same'};
namea={'V1','V4'};
namesign={'minus','plus'};
namep={'tar','test'};

task=['compute LDR for sign ', namea{ba}, ' K = '  sprintf('%1.0i',K)];
display(task)


% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau=0:L;                                                                              % support
lambda=1/tau_prime;                                                                   % time constant
kernel=exp(-lambda.*tau); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

if place ==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results
      
loadname=['weight_bac_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_', sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

%% reconstruct signal in layers

nbses=length(w_all);
ncv=size(w_all{1},1);


x_minus=zeros(nbses,nperm,2,K);
x_plus=zeros(nbses,nperm,2,K);

tic
parfor sess = 1:nbses
    %display(sess)
    %%
    weights=w_all{sess};
    N=size(weights,2);
    
    % get sign
    tell_sign=cell(ncv,2);
    for cv=1:ncv
        tell_sign{cv,1}=find(weights(cv,:)<0);                                          % get the index of positive and negative weights 
        tell_sign{cv,2}=find(weights(cv,:)>0);
    end
    %%
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);             % both conditions
    
    %% randomly permute the label "same" and "different" for spike trains of neurons of the undesired sign
    
    for r = 1:2
        
        %%
        x_perm=zeros(nperm,2,K);
        for perm=1:nperm
            
            rp=randperm(sum(J));                                                                       % random order of trials
            spikes_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false);                     % random order to all spike trains    
            
            for cv=1:ncv
                idx_sign=tell_sign{cv,r};
                spikes_perm{cv}(:,idx_sign,:)=spike_one{cv}(:,idx_sign,:);                         % neurons with sign r (1=neg, 2=pos) have regular spike counts
            end
            
            spikes=cell(ncv,2);
            spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spikes_perm, 'UniformOutput', false);               % put in form cell{ncv,2}
            spikes(:,2)=cellfun(@(x) x(J(1)+1:end,:,:), spikes_perm, 'UniformOutput', false);
            
            %%
            
            [x_cv] = reconstruct_fun(weights,spikes,kernel);
            x_perm(perm,:,:)=x_cv;                                                                      % collect across permutations
            
        end
        
        if r==1
            x_minus(sess,:,:,:)=x_perm;                                                  % average across permutations and collect results across sessions
        else
            x_plus(sess,:,:,:)=x_perm;
        end
        
        
    end
    
end
toc
%%

x_sign=cell(2,1);   % put in a cell
x_sign{1}=x_minus; 
x_sign{2}=x_plus;


%% save

if saveres==1
    savename=['sign_remove_info_',sprintf('%1.0i',T)];
    savefile= '/home/veronika/synced/transfer_result/signal/sign/';
    save([savefile,savename],'x_sign')
    
    clear variables
    
end
        
       
