% reconstruct with random weights an spike timing (jitter in time)

close all
clear all 
clc

place=0;

if place==0
    saveres=1;
else
    saveres=0;
end

nperm=1000;

ba=1;
period=2;
T=400;
 
%%

switcht=139;                                    % flip of the choice signal

nameperms={'entire_window', 'half_window'};
 
perm_type=2;
nameperm=nameperms{perm_type};

namebeh={'same','different'};
namea={'V1','V4'};
namep={'tar','test'};

task=['reconstruct with permuting the spike timing in ',nameperm];
disp(task)

L=100;                              % length of the kernel                                  % exponential kernel for convolution
tau_vec=0:L;                        % support
lambda=1/20;                        % time constant
kernel=exp(-lambda.*tau_vec); 

%% load weights and spikes

addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/') 
end

loadname=['weight_bac_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',T)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_', sprintf('%1.0i',T)];
load(loadname2,'spikes_inm','spikes_cnm')

ncv=size(spikes_cnm,2);
nbses=length(w_all);

%% reconstruct the signal with permutation

tic

K=size(spikes_cnm{1},3);

xp_diff=zeros(nbses,nperm,K);
xp_same=zeros(nbses,nperm,K);
    
for sess=1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spikes=cat(2,spike_diff,spike_same); 
    
    weights=squeeze(w_all{sess});
    
    N=size(spike_diff{1},2);
    J=[size(spike_diff{sess,1},1), size(spike_same{sess,1},1)];
    
    x_perms=zeros(nperm,2,K);
    for perm=1:nperm
           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spikes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
        if perm_type==1 % permute across the entire time window
            for c=1:2
                random_timing=randperm(K);                                                             % permute spike timing vector
                spikes(:,c)=cellfun(@(x) x(:,:,random_timing),spikes(:,c),'UniformOutput',false);
            end
        end
        if perm_type==2 % permute across the entire time window
            for c=1:2
                rt1=randperm(switcht);                                              % permute spike timing vector
                rt2=switcht+randperm(T-switcht);
                rt=cat(2,rt1,rt2);
                spikes(:,c)=cellfun(@(x) x(:,:,rt),spikes(:,c),'UniformOutput',false);
            end
        end
            
        
        %%%%%%%%%%%%%%%%%%%%%% compute LDR %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [x_cv] = reconstruct_fun(weights,spikes,kernel);
        x_perms(perm,:,:) = x_cv;                    
         
    end
    
    xp_diff(sess,:,:)=x_perms(:,1,:);
    xp_same(sess,:,:)=x_perms(:,2,:);
    
end
toc


%% save

if saveres==1
        
    savename=nameperm;
    savefile='/home/veronika/synced/transfer_result/signal/permute_timing/';
    save([savefile,savename],'xp_diff','xp_same')
    %clear all 
    
end
    



