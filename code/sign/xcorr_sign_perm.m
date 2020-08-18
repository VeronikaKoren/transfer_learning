% xcorr sign with permutations

clear all 
close all
clc

place=1;

nperm=1000;

if place==0
    saveres=1;
else
    saveres=0;
end
   
ba=1;
period=2;
tau=20;                                                                % choose between [10,20,50]
T=400;

%%

namea={'V1','V4'};
namep={'tar','test'};

task=['compute xcorr sign permuted ', ' tau='  sprintf('%1.0i',tau)];
display(task)

blue=[0,0.48,0.74];
col={blue,'g'};

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
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

% load results
      
loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',T)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
load(loadname2,'spikes_inm','spikes_cnm')


%% correlation function between signals of plus and minus neurons

nbses=length(w_all);
ncv=size(w_all{1},1);

r_perm=zeros(nbses,nperm,2*T-1);

tic    
parfor sess=1:nbses
    
    wp=w_all{sess};                                                                 % regular weights
    w_abs=max(max(abs(wp)));                                                        % maximal range of weights
    w_pos=wp.*sign(wp);                                                               % makes all weights positive  
    
    N=size(wp,2); 					
    np=floor(N/2);
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spikes_one=cellfun(@(x,y) double(cat(1,x,y)),spike_diff, spike_same,'UniformOutput', false); 
    
    J=size(spikes_one{1},1);                                                            % nb trials
    
    %% compute LDR for plus and minus neurons in trials, with permutation
    
    rxy_perm=zeros(nperm,2*T-1);
    for perm=1:nperm
        
        xp=cell(2,1);                                                                 %(sign, 1)
        wrs=sign(- w_abs + (w_abs + w_abs).*rand(N,100)');                            % random sign of weights
        weights=w_pos.*wrs;                                                           % weights with random sign 
        
        for sgn=1:2
            
            rp=randperm(N);        
            spikes=cellfun(@(x) x(randperm(J),:,:), spikes_one,'UniformOutput',false);
            
            [x_rec] = reconstruct_1c_fun(weights,spikes,kernel);                              % compute reconstruction in trials
            xp{sgn}=x_rec;                                                                    % collect across permutations
           
        end
        
        %% compute correlation function between plus and minus neurons
        
        x=xp{1}; % minus
        y=xp{2}; % plus
        
        [rxy] = correlation_fun(x,y); % correlation between plus and minus signals in trials and cv, averged across trials and cv
        rxy_perm(perm,:)=rxy;
        
    end
    
    r_perm(sess,:,:)=rxy_perm;
    
end
toc

%% save
    
if saveres==1
        
    savename=['xcorr_perm_tau',sprintf('%1.0i',tau),'_',sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/signal/sign/';
    save([savefile,savename],'r_perm')
    clear variables 
end
    


%%

