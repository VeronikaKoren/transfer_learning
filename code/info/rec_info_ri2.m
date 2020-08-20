
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

tau=20;                                                                % choose between [10,20,50]
K=400;

%%
ba=1;
period=2;

namea={'V1','V4'};
namep={'tar','test'};
namebeh={'different','same'};
nameg={'info','noninfo'};

ng=length(nameg);
task=['compute signal ri for neurons with strong and weak weights ', ' K = '  sprintf('%1.0i',K)];
display(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/transfer_mc/'); % for spikes                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/'); % for weights
addpath('/home/veronika/synced/transfer_result/weight/tags/');       % for tag    

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

loadname3=['tag_info_s+c_', sprintf('%1.0i',K)];
load(loadname3)

%% reconstruct signal bursty

nbses=length(w_all);
ncv=size(w_all{1},1);
tag=[1,0];                           % info/noninfo

x_all=zeros(nbses,ng,2,K);           % (sessions, groups, behavior,time)

tic
for sess = 1%:nbses
    
    %%
    weights=w_all{sess};
    N=size(weights,2);
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);             % both conditions
    
    %% randomly permute the label "same" and "different" for spike trains of neurons of other sign
    
    x_partial=zeros(ng,nperm,2,K);
    for r = 1:ng
        
        tagged=find(tag_info{sess}==tag(r));
        [x_cvpb] = reconstruct_partialb_fun(weights,spike_one,kernel,J,nperm,tagged);
        x_partial(r,:,:,:)=x_cvpb;                                                                  
        
    end
    
    x_all(sess,:,:,:)=squeeze(nanmean(x_partial,2)); % average across permutations
    
end
toc


%% save

if saveres==1
    savename=['info_remove_info_',namei{info_case},'_',sprintf('%1.0i',K)];
    savefile= '/home/veronika/synced/transfer_result/signal/info/';
    save([savefile,savename],'x_all')
    
    clear variables
    
end
        
       
