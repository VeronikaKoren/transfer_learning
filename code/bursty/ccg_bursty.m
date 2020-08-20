
% compute ccg noise
% within group
% cross-correlation function (ccg_raw - ccg_trial_invariant)
% 1 condition (concatenated trials from two conditions)

close all
clear all
clc 
format long

saveres=0;

nperm=1000;                                                                       % permutation of sign
nshuffle=20;                                                                    % trial shuffle to subtract the signal correlation

%%

ba=1;
period=2;

K=400;
start_vec=[200,500];                                                             % beginning of the time window for the target (200) and the test stimulus (500) 
start=start_vec(period);
taumax=50;

namea={'V1','V4'};
namep={'target', 'test'};
nameg={'bursty','nonbursty'};

%% load spike counts                                                                         

addpath('/home/veronika/synced/transfer_result/weight/tags/');
addpath('/home/veronika/synced/transfer_result/input/spike_train/');
addpath('/home/veronika/Dropbox/transfer/code/function/')

loadname=['spike_train_c_',namea{ba},'_',namep{period}];   % spike trains from choice (CNM/INM)
load(loadname);

strain=cellfun(@(x,y) single(cat(1,x(:,:,start:start+K-1),y(:,:,start:start+K-1))),spiketrain(:,1),spiketrain(:,2), 'UniformOutput', false);

%% load weights

loadname3=['tag_bursty_', sprintf('%1.0i',K)]; % tag bursty from S+C
load(loadname3,'bursty')

nbses=length(strain);
%%                                                                          

g1=cell(nbses,1);
g2=cell(nbses,1);

gp1=cell(nbses,nperm);
gp2=cell(nbses,nperm);

tag=[1,0];

tic
parfor sess = 1:nbses
    
    st=strain{sess};
    N=size(st,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regular
    % bursty/bursty
    idx1=find(bursty{sess}==tag(1));
    if length(idx1)>1
        spike_train=st(:,idx1,:);
        g1{sess} = ccg_fun(spike_train, nshuffle);
    end
    
    % nonbursty/nonbursty
    idx2=find(bursty{sess}==tag(2));
    if length(idx2)>1
        spike_train=st(:,idx2,:);
        g2{sess} = ccg_fun(spike_train, nshuffle);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% permuted
    %%
    
    for perm=1:nperm
        rndb=binornd(1,0.5,N,1);                     % random draw of N [0,1] from binomial distribution
        idx1=find(rndb==tag(1));
        if length(idx1)>1
            spike_train=st(:,idx1,:);
             ccg= ccg_fun(spike_train, nshuffle);
             gp1{sess,perm}=ccg(:,K-taumax:K+taumax);
        end
        
        % nonbursty/nonbursty
        idx2=find(rndb==tag(2));
        if length(idx2)>1
            spike_train=st(:,idx2,:);
             ccg= ccg_fun(spike_train, nshuffle);
             gp2{sess,perm}=ccg(:,K-taumax:K+taumax);
        end
        
    end
    
end
toc

%%
lags=-K+1:K-1;

g1mat=cell2mat(g1);
g2mat=cell2mat(g2);

%% permutation test for difference in synchrony between plus and minus

d=mean(g1mat(:,K))-mean(g2mat(:,K));

% sum pairs across sessions
p1=cell(nperm,1);
p2=cell(nperm,1);
for perm=1:nperm
    p1{perm}=cell2mat(gp2(:,perm));              
    p2{perm}=cell2mat(gp1(:,perm));
end

miniK=taumax;
d0=cellfun(@(x,y) mean(x(:,miniK))-mean(y(:,miniK)),p1,p2); % compute the synchrony, averaged across pairs
pval=sum(d<d0)/nperm;
display(pval,'permutation test group 1 stronger than group 2')

%% permutation test for rccg for different summation length

ntau=20;

dc=zeros(ntau,1);
d0mm=zeros(ntau,2);
pval_coeff=zeros(ntau,1);

r1=zeros(ntau,size(g1mat,1));
r2=zeros(ntau,size(g2mat,1));

for tau=1:ntau
    
    rc1=sum(g1mat(:,K-tau:K+tau),2);
    rc2=sum(g2mat(:,K-tau:K+tau),2);
    
    r1(tau,:)=rc1;
    r2(tau,:)=rc2;
    
    dc(tau)=mean(rc1) - mean(rc2);
    dc0=cellfun(@(x,y) mean(sum(x(:,miniK-tau:miniK+tau),2))-mean(sum(y(:,miniK-tau:miniK+tau),2)),p1,p2);
    pval_coeff(tau)=sum(dc(tau)<dc0)/nperm;
    
    d0mm(tau,1)=min(dc0);
    d0mm(tau,2)=max(dc0);
    
end

%% save results

if saveres==1
    address='/home/veronika/synced/transfer_result/pairwise/ccg/';
    filename=['ccg_bursty_',sprintf('%1.0i',K)];
    save([address, filename], 'lags','g1','g2','r1','r2','gp1','gp2','pval','pval_coeff','tauvec')
    %clear all
end





