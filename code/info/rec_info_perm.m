% reconstruct layers

clear all 
close 
clc
format short

place=1;

if place==0
    saveres=1;
else
    saveres=1;
end
 
nperm=1000;

ba=1;
period=2;
tau=20;                                                                 
K=400;
%%

namea={'V1','V4'};
namep={'tar','test'};

task='compute signal info permuted ';
disp(task)

% exponential kernel for convolution
L=100;                                                                       % length of the kernel                                  
tau_vec=0:L;                                                                 % support
lambda=1/tau;                                                                % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end
addpath('/home/veronika/synced/transfer_result/input/transfer_mc/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');
addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/');
addpath('/home/veronika/synced/transfer_result/weight/tags/');

%% load results

loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')
      
load('tag_info_s+c_400')

nbses=length(w_all);
ncv=size(w_all{1},1);

%%

xp1=zeros(nbses,nperm,2,K);
xp2=zeros(nbses,nperm,2,K);


tic
for sess = 1:nbses
    %disp(sess)
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);             % concatenate conditions
    
    J=[size(spike_diff{sess,1},1), size(spike_same{sess,1},1)];
    N=size(spike_diff{1},2);
    %%
    w_regular=squeeze(w_all{sess});
    w_abs=max(max(abs(w_regular)));
    
    x_1=zeros(nperm,2,K);
    x_2=zeros(nperm,2,K);
    
    for perm=1:nperm                                                                             % permutation cycles
        
        w_random=-w_abs + (w_abs + w_abs).*rand(1,N);
        
        spike_perm=cellfun(@(x) x(randperm(sum(J)),:,:),spike_one, 'UniformOutput',false);        % permute class label for spike trains
        spikes=cell(ncv,2);
        spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spike_perm, 'UniformOutput',false);
        spikes(:,2)=cellfun(@(x) x(J(1)+1:sum(J),:,:), spike_perm, 'UniformOutput',false);
    
        for r = 1:2
                
            weights=repmat(w_random,ncv,1);
            [x_cv] = reconstruct_fun(weights,spikes,kernel);
            
            if r==1
                x_1(perm,:,:)=x_cv;                                                 % collect results across permutation
            elseif r==2
                x_2(perm,:,:)=x_cv;
            end
            
        end
    end
    
    xp1(sess,:,:,:)=x_1;                                                % collect results across sessions
    xp2(sess,:,:,:)=x_2;
    
end

toc

%%

xperm=cell(2,1);   % put in a cell
xperm{1}=xp1; 
xperm{2}=xp2;

%% save

if saveres==1
    
    savename='info_perm';
    savefile= '/home/veronika/synced/transfer_result/signal/info/';
    save([savefile,savename],'xperm')
    clear all
    
end
        

%exit
