
% compute ccg noise within layers 

close all
clear all
clc 
format long

saveres=0;
info_case=2;                                                                    % determines the info case for weights (use S+C or C)

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
namei={'s+c','c'};
nameg={'SG','G','IG'};
ng=length(nameg);


%% load                                                                          

addpath('/home/veronika/synced/transfer_result/weight/tags/');
addpath('/home/veronika/synced/transfer_result/input/spike_train/');
addpath('/home/veronika/Dropbox/transfer/code/function/')
addpath('/home/veronika/synced/transfer_result/signal/layer/');

loadname=['spike_train_',namei{info_case},'_',namea{ba},'_',namep{period}];                            % spike trains 
load(loadname);

strain=cellfun(@(x,y) single(cat(1,x(:,:,start:start+K-1),y(:,:,start:start+K-1))),spiketrain(:,1),spiketrain(:,2), 'UniformOutput', false);

loadname3=['ncell_lay_',namea{ba}];
load(loadname3)

nbses=length(strain);

%%                                                                          

ccg3=cell(3,1);
perm3=cell(3,1);

tic

for g=1:3
    
    ccg_sess=cell(nbses,1);
    ccgp_sess=cell(nbses,nperm);
    
    parfor sess = 1:nbses
        
        st=strain{sess};
        N=size(st,2);
        
        ncl=ncell_lay(sess,:);
        s=cumsum([0,ncl]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% regular
        
        tagged = s(g) + 1 : s(g+1);
        spike_train=st(:,tagged,:);
        ccg_sess{sess}=ccg_fun(spike_train, nshuffle);
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% permuted
          
        for perm=1:nperm
            
            rndb=randi(3,[N,1])
            idx=find(rndb==g);
            
            if length(idx)>1
                
                spike_train=st(:,idx,:);
                ccg= ccg_fun(spike_train, nshuffle);
                ccgp_sess{sess,perm}=ccg(:,K-taumax:K+taumax);
      
            end
        end  
    end
    %%
    ccg3{g}=cell2mat(ccg_sess);
    
    ccgp=cell(nperm,1);
    for perm=1:nperm
        ccgp{perm}=cell2mat(ccgp_sess(:,perm));    % collect across sessions, for every permutation cycle
    end
    
    perm3{g}=ccgp;
    
end
toc

lags=-K+1:K-1;

%% permutation test for difference in synchrony across layers

syn=cellfun(@(x) x(:,K),ccg3,'UniformOutput', false);       % synchrnony
msyn=cellfun(@mean, syn);                                   % average across pairs

idx1=[1,1,2];
idx2=[2,3,3];

pval=zeros(3,1);
for ii=1:3
    
    d=msyn(idx1(ii))-msyn(idx2(ii));
    
    p1=perm3{idx1(ii)};
    p2=perm3{idx2(ii)};
    
    synp1=cellfun(@(x) mean(x(:,taumax)), p1);
    synp2=cellfun(@(x) mean(x(:,taumax)), p2);
    
    d0=synp1-synp2;
    
    pval(ii)=sum(d<d0)/nperm;
   
end

%% permutation test rccg for different summation length

ntau=20;
r_all=cell(3,1);
for g=1:3
    
    x=ccg3{g};
    r=zeros(ntau,size(x,1));
    for tau=1:ntau
        r(tau,:)=sum(x(:,K-tau:K+tau),2);
    end
    r_all{g}=r;
end
%%

pval_coeff=zeros(3,ntau);

for ii=1:3
    
    x=r_all{idx1(ii)};
    y=r_all{idx2(ii)};
    
    p1=perm3{idx1(ii)};
    p2=perm3{idx2(ii)};
    
    
    for tau=1:ntau
        
        dc=mean(x(tau,:)) - mean(y(tau,:));
        
        synp1=cellfun(@(x) mean(sum(x(:,miniK-tau:miniK+tau),2)), p1);
        synp2=cellfun(@(x) mean(sum(x(:,miniK-tau:miniK+tau),2)), p2);
        
        dc0=synp1-synp2;
        pval_coeff(ii,tau)=sum(dc < dc0)/nperm;
        
    end
end
%%

tauvec=1:ntau;
idxes=[{idx1},{idx2}];

%% save results

if saveres==1
    address='/home/veronika/synced/transfer_result/pairwise/ccg/';
    filename=['ccg_layer_',namei{info_case},'_',sprintf('%1.0i',K)];
    save([address, filename], 'lags','ccg3','perm3','r_all','pval','pval_coeff','tauvec','idxes')
    %clear all
end

