% correlation function between plus and minus neurons

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

nperm=1000;       
ba=1;
period=2;

tau=20;                                                               
K=400;

%%

namea={'V1','V4'};
namep={'tar','test'};
namelay={'SG','G','IG'};

name2lay=cell(3,1);
idx1=[1,1,2];
idx2=[2,3,3];
for i=1:3
    
    lay1=namelay{idx1(i)};
    lay2=namelay{idx2(i)};
    
    name2lay{i}=[lay1,' & ',lay2];
end

task=['compute xcorr for pairs of layers ', ' K='  sprintf('%1.0i',K)];
display(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                          % support
lambda=1/tau;                                                                         % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/reconstruction/result/layer/ncell/')
if place==1
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

%% correlation function between signals of plus and minus neurons

nbses=length(w_all);
ncv=size(w_all{1},1);
lags=-K+1:K-1;

tic
r_lay_perm=cell(nbses,1);

parfor sess=1:nbses
   
    w_sess=w_all{sess};
    w_pos=w_sess.*sign(w_sess);
    
    ncl=ncell_lay(sess,:);
    s=cumsum([0,ncl]);
    N=size(w_sess,2);
    
    w_abs=max(max(abs(w_sess)));
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) double(cat(1,x,y)),spike_diff, spike_same,'UniformOutput', false);                                                
    J=size(spike_one{1},1);
    
    %% compute xcorr layer perm
    
    r_perm=zeros(3,nperm,2*K-1);
    for perm=1:nperm                                                                             % permutation cycles
        
        
        xp=cell(3,1);
        
        for lay = 1:3
            
            w_random=sign(-w_abs + (w_abs + w_abs).*rand(1,N));
            weights=w_random.*w_pos;
            
            spikes=cellfun(@(x) x(randperm(J),:,:),spike_one, 'UniformOutput',false);
            [x_rec] = reconstruct_1c_fun(weights,spikes,kernel);                               % compute reconstruction in trials
            
            xp{lay}=x_rec;                                                                   % collect across permutations
            
        end
        
        %% compute correlation function between the population signal in layer 1 and layer 2
        for i=1:3
            
            x=xp{idx1(i)}; % layer 1
            y=xp{idx2(i)}; % layer 2
            
            [rxy] = correlation_fun(x,y); % correlation between pop signal for a pair oflayers, averged across trials and cv
            r_perm(i,perm,:)=rxy;
            
        end
        
    end
    
    r_lay_perm{sess}=r_perm;
    
end
toc

%% convert into a matrix and average across sessions
if nbses>2
    r_layer_perm=zeros(nbses,3,nperm,2*K-1);
    for sess=1:nbses
        r_layer_perm(sess,:,:,:)=r_lay_perm{sess};
    end
end
%% plotfig

if pltfig==1
    rmean=squeeze(nanmean(r_layer_perm));
    figure()
    hold on
    for perm=1:nperm
        
        plot(squeeze(rmean(1,perm,:)),'k')
        
    end
end
    
%% save

if saveres==1
        
    savename=['xcorr_layer_perm_tau',sprintf('%1.0i',tau)];
    savefile='/home/veronika/synced/transfer_result/signal/layer/';
    save([savefile,savename],'r_layer_perm','lags','name2lay')
    clear variables
    
end
    