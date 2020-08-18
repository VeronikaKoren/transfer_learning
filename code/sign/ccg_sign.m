
% compute ccg noise
% within group (within plus and within minus neurons) 
% cross-correlation function (ccg_raw - ccg_trial_invariant)
% 1 condition (concatenated trials from two conditions)

close all
clear all
clc 
format long

saveres=0;
showfig=1;
info_case=2;                                                                    % determines the info case for weights (use S+C or C)

nperm=1000;                                                                       % permutation of sign
nshuffle=20;                                                                    % trial shuffle to subtract the signal correlation

%%

ba=1;
period=2;

K=400;
start_vec=[200,500];                                                             % beginning of the time window for the target (200) and the test stimulus (500) 
start=start_vec(period);

namea={'V1','V4'};
namep={'target', 'test'};
namei={'s+c','c'};

%% load spike counts                                                                         

addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/');
addpath('/home/veronika/synced/transfer_result/input/spike_train/');
addpath('/home/veronika/Dropbox/transfer/code/function/')

loadname=['spike_train_',namea{ba},'_',namep{period}];                            % spike trains from choice (CNM/INM)
load(loadname);

strain=cellfun(@(x,y) single(cat(1,x(:,:,start:start+K-1),y(:,:,start:start+K-1))),spiketrain(:,1),spiketrain(:,2), 'UniformOutput', false);

%% load weights

loadname2=['weight_alltr_',namei{info_case},'_',namea{ba},'_',sprintf('%1.0i',K),'.mat'];    % use weights from s+c                
load(loadname2)

nbses=size(w_alltr,1);
%%                                                                          

cminus=cell(nbses,1);
cplus=cell(nbses,1);

pp=cell(nbses,nperm);
pm=cell(nbses,nperm);

tic
parfor sess = 1:nbses
    
    st=double(strain{sess});
    w=w_alltr{sess};
    N=length(w);
    w_abs=max(abs(w));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regular
    % (-/-)
    neg=find(w<0);
    if length(neg)>1
        spike_train=st(:,neg,:);
        cminus{sess} = ccg_fun(spike_train, nshuffle);
    end
    
    % (+/+)
    pos=find(w>0);
    if length(pos)>1
        spike_train=st(:,pos,:);
        cplus{sess} = ccg_fun(spike_train, nshuffle);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% permuted
    %%
    
    w_nosign=w.*sign(w);                                                    % all positive
    
    for perm=1:nperm
        wrs=sign(- w_abs + (w_abs + w_abs).*rand(N,1)')';                    % random sign
        w_perm=w_nosign.*wrs;                                                % apply random sign
        
        % p1
        neg=find(w_perm<0);
        
        if length(neg)>1
            spike_train=st(:,neg,:);
            [ccg]=ccg_fun(spike_train, nshuffle);
            pm{sess,perm}=ccg(:,K-50:K+50);
        end
        
        % p2
        pos=find(w_perm>0);
        
        if length(pos)>1
            spike_train=st(:,pos,:);
            [ccg]=ccg_fun(spike_train, nshuffle);
            pp{sess,perm}=ccg(:,K-50:K+50);
        end
       
    end
   
end
toc

%%
lags=-K+1:K-1;

ccgplus=cell2mat(cplus);
ccgminus=cell2mat(cminus);

%% permutation test for difference in synchrony between plus and minus

d=mean(ccgplus(:,K))-mean(ccgminus(:,K));

% sum pairs across sessions
p1=cell(nperm,1);
p2=cell(nperm,1);
for perm=1:nperm
    p1{perm}=cell2mat(pm(:,perm));              
    p2{perm}=cell2mat(pp(:,perm));
end

miniK=(size(p1{1},2)-1)/2;
d0=cellfun(@(x,y) mean(x(:,miniK))-mean(y(:,miniK)),p1,p2); % compute the synchrony, averaged across pairs
pval=sum(d<d0)/nperm;
display(pval,'permutation test')

%% permutation test for rccg for different summation length

ntau=20;

dc=zeros(ntau,1);
d0mm=zeros(ntau,2);
pval_coeff=zeros(ntau,1);

r1=zeros(ntau,size(ccgplus,1));
r2=zeros(ntau,size(ccgminus,1));

for tau=1:ntau
    
    rc1=sum(ccgplus(:,K-tau:K+tau),2);
    rc2=sum(ccgminus(:,K-tau:K+tau),2);
    
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
    filename=['ccg_sign_',namei{info_case},'_',namea{ba},namep{period}];
    save([address, filename], 'lags','cplus','cminus','r1','r2','pp','pm','pval','pval_coeff','tauvec')
    %clear all
end





