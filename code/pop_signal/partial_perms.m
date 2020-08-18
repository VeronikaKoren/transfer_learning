% reconstruct with random weights and with permuted class for spikes

close all
clear all 
clc

place=1;

if place==0
    showfig=0;
    saveres=1;
else
    saveres=0;
    showfig=1;
    factor=100;
end

nperm=1000;

ba=1;
period=2;
T=400;
 
%%

choose_perm={'1:random_wsign','2:random_wamp','3:permute_nspikes','4:permute_cspikes',};
nameperms=cellfun(@(x) x(3:end),choose_perm,'UniformOutput',false);

perm_type=1;
nameperm=nameperms{perm_type};

namebeh={'same','different'};
namea={'V1','V4'};
namep={'tar','test'};

task=['reconstruct with ',nameperm];
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
    
parfor sess=1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    
    if perm_type==4     % permute class for spikes
        spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);
    else
        spikes=cat(2,spike_diff,spike_same);
    end
    
    w_regular=squeeze(w_all{sess});
    w_abs=max(max(abs(w_regular)));
    
    N=size(spike_diff{1},2);
    J=[size(spike_diff{sess,1},1), size(spike_same{sess,1},1)];
    
    x_perms=zeros(nperm,2,K);
    for perm=1:nperm
           
        %%%%%%%%%%%%%%%%%%%%%%%% weights
        if perm_type>2
            weights=w_regular;
        else
            w_random=repmat(-w_abs + (w_abs + w_abs).*rand(1,N),ncv,1);
            
            if perm_type==1
                weights=(w_regular.*sign(w_regular)).*sign(w_random);       % random sign, regular amplitude 
            elseif perm_type==2
                weights=(w_random.*sign(w_random)).*sign(w_regular);        % random amplitude, regular sign
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spikes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%
        
        if perm_type==3 % permute spikes across neurons
            rp=randperm(N);
            spikes=cellfun(@(x) x(:,rp,:), spikes, 'UniformOutput', false);
            
        end
        
        if perm_type==4 % permute class for spikes
            spike_perm=cellfun(@(x) x(randperm(sum(J)),:,:),spike_one, 'UniformOutput',false);        % permute class label for spike trains
            spikes=cell(ncv,2);
            spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spike_perm, 'UniformOutput',false);
            spikes(:,2)=cellfun(@(x) x(J(1)+1:sum(J),:,:), spike_perm, 'UniformOutput',false);
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
    savefile='/home/veronika/synced/transfer_result/signal/permuted/';
    save([savefile,savename],'xp_diff','xp_same')
    %clear all 
    
end
    

