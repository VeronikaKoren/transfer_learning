% reconstruct layers

clear all 
close 
clc
format short

place=1;

if place==0
    saveres=1;
else
    saveres=0;
end
 
nperm=1000;

ba=1;
period=2;

tau=20;                                                                 
K=400;
%%

namea={'V1','V4'};
namelay={'SG','G','IG'};
namep={'tar','test'};

task='compute signal in layers permuted ';
disp(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end
addpath('/home/veronika/synced/transfer_result/input/transfer_mc/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/synced/transfer_result/signal/layer/');

%% load results
      
loadname=['weight_bac_',namea{ba},'_',namep{period}, '_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_' sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['ncell_lay_',namea{ba}];
load(loadname3);

nbses=length(w_all);
ncv=size(w_all{1},1);

%%

x_sgp=zeros(nbses,nperm,2,K);
x_gp=zeros(nbses,nperm,2,K);
x_igp=zeros(nbses,nperm,2,K);

tic
for sess = 1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    J=[size(spike_diff{sess,1},1), size(spike_same{sess,1},1)];
    
    ncl=ncell_lay(sess,:);
    s=cumsum([0,ncl]);
    N=s(end);
    
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);             % concatenate conditions
    
    %%
    w_regular=squeeze(w_all{sess});
    w_abs=max(max(abs(w_regular)));
    
    x_sg=zeros(nperm,2,K);
    x_g=zeros(nperm,2,K);
    x_ig=zeros(nperm,2,K);
    
    for perm=1:nperm                                                                             % permutation cycles
        
        w_random=-w_abs + (w_abs + w_abs).*rand(1,N);                                            % random weights   
        
        spike_perm=cellfun(@(x) x(randperm(sum(J)),:,:),spike_one, 'UniformOutput',false);       % permute class label for spike trains
        spikes=cell(ncv,2);
        spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spike_perm, 'UniformOutput',false);
        spikes(:,2)=cellfun(@(x) x(J(1)+1:sum(J),:,:), spike_perm, 'UniformOutput',false);
    
        for r = 1:3
                
            weights=repmat(w_random,ncv,1);
            [x_cv] = reconstruct_fun(weights,spikes,kernel);
            
            if r==1
                x_sg(perm,:,:)=x_cv;                                                 % collect results across permutation
            elseif r==2
                x_g(perm,:,:)=x_cv;
            elseif r==3
                x_ig(perm,:,:)=x_cv;
            end
            
        end
        
    end
    % perm
    
    x_sgp(sess,:,:,:)=x_sg;                                                % collect results across sessions
    x_gp(sess,:,:,:)=x_g;
    x_igp(sess,:,:,:)=x_ig;
    
end

toc
%%
xp_layer=cell(3,1);
xp_layer{1}=x_sgp;
xp_layer{2}=x_gp;
xp_layer{3}=x_igp;

%% save

if saveres==1
    
    savename='layer_permuted';
    savefile= '/home/veronika/synced/transfer_result/signal/layer/';
    save([savefile,savename],'xp_layer')
    %clear all
    
end
        
disp('finished')

