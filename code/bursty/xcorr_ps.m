
% correlation function between info and noninfo neurons
% for all trials

clear all 
close all
clc

place=1;

if place==1
    saveres=0;
    pltfig=1;
else
    saveres=1;
    pltfig=0;
end

nperm=3;

ba=1;
period=2;
tau=20;                                                                % choose between [10,20,40]
T=400;

%%

namea={'V1','V4'};
namep={'tar', 'test'};

task=['compute xcorr between bursty and non-bursty, 1c ', ' T='  sprintf('%1.0i',T), ' tau ='  sprintf('%1.0i',tau)];
display(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/synced/transfer_result/ps/');

if place==1
    addpath('/home/veronika/Dropbox/matlabfun/shadedErrorBar')
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results
      
loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',T)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['tag_bursty_',namep{period},namea{ba}];
load(loadname3,'bursty')

%% compute

nbses=length(w_all);
ncv=size(w_all{1},1);
lags=-T+1:T-1;

r_ps=zeros(nbses,2*T-1);

tic
parfor sess=1:nbses
   
    weights=w_all{sess};
    
    %% get tags
    tag=bursty{sess};                         % long ISI
    anti_tag=(tag==0);                      % short ISI
    tags=cat(2,tag,anti_tag);
    
    %%
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y),spike_diff, spike_same,'UniformOutput', false);                                                
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    %% compute LDR in trials
    
    r_perm=zeros(nperm,2*T-1);
    for perm=1:nperm
        
        xp=cell(2,1);
        for r=1:2
            
            use_tags=double(find(tags(:,r)));
            
            rp=randperm(sum(J));                                                                       % random order of trials
            spikes_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false);                     % random order to all spike trains    
            
            for cv=1:ncv
                spikes_perm{cv}(:,use_tags,:)=spike_one{cv}(:,use_tags,:);                             % tagged neurons have correct label
            end
            
            spikes=spikes_perm;
            
            %%
            
            [x_rec] = reconstruct_1c_fun(weights,spikes,kernel);
            xp{r}=x_rec;                                                                             
                                                                        
        end
        
        %% compute correlation function 
    
        x=xp{1}; % group 1
        y=xp{2}; % group 2
    
        [rxy] = correlation_fun(x,y); % correlation between population signals of group 1 and group 2 in trials and cv, averged across trials and cv
        r_perm(perm,:)=rxy;
        
    end
    
    r_ps(sess,:)=nanmean(r_perm);
    
end
toc

%%
if pltfig==1
    
    figure()
    hold on
    for sess=1:nbses
        subplot(4,5,sess)
        plot(r_ps(sess,:))
        
    end
    
    figure()
    plot(mean(r_ps))
  
end

%%

%% save

if saveres==1
        
    savename=['xcorr_ps_tau',sprintf('%1.0i',tau)];
    savefile='/home/veronika/synced/transfer_result/signal/ps/';
    save([savefile,savename],'r_isi','lags')
    %clear variables
    
end
    
%exit

%%

