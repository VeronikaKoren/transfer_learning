

close all
clear all
clc

saveres=0;

%%
ba=1;
period=2;

namep={'tar','test'};
namea={'V1','V4'};
namei={'s+c','c'};

K=400;          % length of the time window

%% load data 

addpath('/home/veronika/synced/transfer_result/basic_stat/ps/')
loadname=['ps_',namea{ba},namep{period}, '_',sprintf('%1.0i',K)];
load(loadname)

loadname2=['ps_poisson_',namea{ba},namep{period}, '_',sprintf('%1.0i',K)];
load(loadname2)

nf=length(f);   % number of frequencies evauated
nbses=size(psp,1);
nperm=size(psp{1},1);

%% tag bursty neurons

vec=3:41; % use frequencies from 10 to 200 Hz;

aps=cellfun(@(x) nanmean(x(:,vec),2), ps_sess,'UniformOutput', false);          % area under the power spectum regular model
aps_perm=cellfun(@(x) nanmean(x(:,:,vec),3), psp,'UniformOutput', false);       % area poisson model

alpha=0.05;
idx=nperm*alpha;

lower_bound=cell(nbses,1);
bursty=cell(nbses,1);
for sess=1:nbses
    
    N=size(ps_sess{sess},1);
    burst=zeros(N,1);
    
    lb=zeros(N,1);
    for i=1:N
       x=aps{sess}(i);
       x0=aps_perm{sess}(:,i);
       
       sorted=sort(x0);
       lb(i)=sorted(idx);
       ub=sorted(end-idx);
       
       burst(i)=x<lb(i);                  % tag=1 for bursty neurons (area under ps is smaller than for poisson)
       
    end
    
    lower_bound{sess}=lb;
    bursty{sess}=burst;
    
end

nbursty=sum(cell2mat(bursty));
display(nbursty,'number of bursty neurons')
%%
used_freq=f(vec);
%% save

if saveres==1

    savename=['tag_bursty_', sprintf('%1.0i',K)];
    savefile='/home/veronika/synced/transfer_result/weight/tags/';
    save([savefile,savename],'bursty','alpha','used_freq','aps','lower_bound')
end


