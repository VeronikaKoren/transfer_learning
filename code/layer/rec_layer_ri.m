% reconstruct layers; remove the information from two out of three layers
% by using random permutation of the label on the spike train. This creates
% a "zero" signal for the two layers and shows the contribution of the
% desired layer.

clear all 
close all
clc

place=1;

if place==1
    saveres=0;
    pltfig=1;
else
    pltfig=0;
    saveres=1;
end

ba=1;
period=2;
tau=20;                                                                % choose between [10,20,40]

K=400;
nperm=1000;
%%

namebeh={'same','different'};
namea={'V1','V4'};
namep={'tar','test'};
namelay={'SG','G','IG'};

task=['compute LDR in layers ', namea{ba}, ' T='  sprintf('%1.0i',K)];
display(task)

blue=[0,0.48,0.74];
col={'g',blue};

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/transfer_mc/')
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/synced/transfer_result/signal/layer/');

if place ==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results
      
loadname=['weight_bac_',namea{ba},'_',namep{period}, '_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_' sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['ncell_lay_',namea{ba}];
load(loadname3);

%% reconstruct signal in layers

nbses=length(w_all);
ncv=size(w_all{1},1);

x_sg=zeros(nbses,nperm,2,K);
x_g=zeros(nbses,nperm,2,K);
x_ig=zeros(nbses,nperm,2,K);


tic
parfor sess = 1:nbses
    
    weights=w_all{sess};
    ncl=ncell_lay(sess,:);
    s=cumsum([0,ncl]);
    N=size(weights,2);
    
    %%
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false); % both conditions
    
    %% randomly permute the label "same" and "different" for spike trains of neurons that are not of interest
    
    x_perm=zeros(nperm,2,K);
    for r = 1:3
        
        delta = s(r) + 1 : s(r+1);
        
        %%
        for perm=1:nperm
            
            rp=randperm(sum(J));                                                                       % random order of trials
            spikes_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false);                     % random order to all spike trains    
            
            
            for cv=1:ncv
                spikes_perm{cv}(:,delta,:)=spike_one{cv}(:,delta,:);                               % neurons from the layer r have regular spike trains
            end
            
            spikes=cell(ncv,2);
            spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spikes_perm, 'UniformOutput', false);               % put in form cell{ncv,2}
            spikes(:,2)=cellfun(@(x) x(J(1)+1:end,:,:), spikes_perm, 'UniformOutput', false);
            
            %%
            
            [x_cv] = reconstruct_fun(weights,spikes,kernel);
            
            x_perm(perm,:,:)=x_cv;                                                                     % collect across perm

        end
        
        if r==1
            x_sg(sess,:,:,:)=x_perm;                                                  % collect results across sessions and average across perm
        elseif r==2
            x_g(sess,:,:,:)=x_perm;
        elseif r==3
            x_ig(sess,:,:,:)=x_perm;
        end
        
    end
    
end
toc

%%

x_lay=cell(3,1);
x_lay{1}=x_sg; 
x_lay{2}=x_g;
x_lay{3}=x_ig;


%% save

if saveres==1
    
    savename='layer_remove_info';
    savefile='/home/veronika/synced/transfer_result/signal/layer/';
    save([savefile,savename],'x_lay')
    
end
%}        
       
