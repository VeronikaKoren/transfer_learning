% reconstruct with random weights and with permuted class for spikes

close all
clear all 
clc

place=1;

if place==0
    saveres=1;
else
    saveres=0;
end

nperm=1000;
info_case=1;  								          % 1 is for transfer of learning, 2 is for choice only

ba=1;
period=2;

T=400;
 
%%

namebeh={'same','different'};
namea={'V1','V4'};
namep={'tar','test'};
namec={'','c'};

nameperm=['permute_class_', namec{info_case}];

task=['reconstruct ',nameperm,'_', namea{ba},'_',namep{period}];
disp(task)

L=100;                          % length of the kernel                                  % exponential kernel for convolution
tau_vec=0:L;                    % support
lambda=1/20;                    % time constant
kernel=exp(-lambda.*tau_vec); 

%% load weights and spikes

addpath('/home/veronika/synced/transfer_result/input/');                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/') 
end

if info_case==1
    loadname=['weight_bac_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',T)];
    load(loadname,'w_all')
    
    loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
    load(loadname2,'spikes_inm','spikes_cnm')
    
else
    
    loadname=['weight_choice_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',T)];
    load(loadname,'w_all')
    
    loadname2=['input_c',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
    load(loadname2,'spikes_inm','spikes_cnm')
    
end


ncv=size(spikes_cnm,2);
nbses=length(w_all);

%% reconstruct the signal with permutation

tic

T=size(spikes_cnm{1},3);

xp_diff=zeros(nbses,nperm,T);
xp_same=zeros(nbses,nperm,T);
    
parfor sess=1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);
    
    w_regular=squeeze(w_all{sess});
    w_abs=max(max(abs(w_regular)));
    
    N=size(spike_diff{1},2);
    J=[size(spike_diff{sess,1},1), size(spike_same{sess,1},1)];
    
    x_perms=zeros(nperm,2,T);
    for perm=1:nperm
       
        weights=repmat(-w_abs + (w_abs + w_abs).*rand(1,N),ncv,1);
        
        spike_perm=cellfun(@(x) x(randperm(sum(J)),:,:),spike_one, 'UniformOutput',false);        % permute class label for spike trains
        spikes=cell(ncv,2);
        spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spike_perm, 'UniformOutput',false);
        spikes(:,2)=cellfun(@(x) x(J(1)+1:sum(J),:,:), spike_perm, 'UniformOutput',false);
        
        [x_cv] = reconstruct_fun(weights,spikes,kernel);
        x_perms(perm,:,:) = x_cv;                    
         
    end
    
    xp_diff(sess,:,:)=x_perms(:,1,:);
    xp_same(sess,:,:)=x_perms(:,2,:);
    
end
toc

%% save result
 
if saveres==1
        
    savename=[nameperm, '_', namea{ba},namep{period}];
    savefile='/home/veronika/synced/transfer_result/signal/permuted/';
    save([savefile,savename],'xp_diff','xp_same')
    clear all 
    
end
    



