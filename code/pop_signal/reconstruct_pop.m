% reconstruct regular with different time constants

clear all 
close all
clc

ba=1;
period=2;

saveres=0;

info_case=1;

tau_range=[10,20,40];
K_range=[300,400,500];

tau=tau_range(2);                                                       %
K=K_range(2);
%%

namebeh={'same','different'};
namea={'V1','V4'};

namep={'tar','test'};
namec={'','c'};

task=['compute LDR all neurons ',' ', namec{info_case},' ',namep{period},' T='  sprintf('%1.0i',K)];
display(task)
                                                                                      % exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                          % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/Dropbox/transfer/code/function/')
addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

%% 
if info_case==1
    loadname=['weight_bac_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',K)];
    load(loadname,'w_all')
    
    loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',K)];
    load(loadname2,'spikes_inm','spikes_cnm')
    
else
    
    loadname=['weight_choice_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',K)];
    load(loadname,'w_all')
    
    loadname2=['input_c',namea{ba},'_', namep{period},'_',sprintf('%1.0i',K)];
    load(loadname2,'spikes_inm','spikes_cnm')
    
end
ncv=size(spikes_cnm,2);
nbses=length(w_all);


%% reconstruct the signal as the weighted sum of spikes

x_diff=zeros(nbses,K);
x_same=zeros(nbses,K);


for sess=1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spikes=cat(2,spike_diff,spike_same);
    
    weights=w_all{sess};
    [x_cv] = reconstruct_fun(weights,spikes,kernel);                        % to reconstruct with regular weights
    
    x_diff(sess,:)=x_cv(1,:);                                                   % collect results across sessions; non-match
    x_same(sess,:)=x_cv(2,:);                                                   % match
    
end

%% save

if saveres==1
    savename=['signal_',namec{info_case},namea{ba},namep{period} sprintf('%1.0i',K),'_tau', sprintf('%1.0i',tau)];
    savefile='/home/veronika/synced/transfer_result/signal/regular/';
    save([savefile,savename],'x_same','x_diff')
    
end
        
%%


