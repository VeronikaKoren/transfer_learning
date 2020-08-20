

clear all 
close all
clc

saveres=0;
ba=1;
T=400;

%%%

namea={'V1','V4'};

task='determine neurons with strong and weak weights';
disp(task)
%%

addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/');

loadname=['weight_alltr_s+c_',namea{ba},'_',sprintf('%1.0i',T)];
load(loadname)

loadname2=['weight_allp_s+c_',namea{ba},'_',sprintf('%1.0i',T)];
load(loadname2)
nbses=size(w_alltr,1);

%%
nperm=size(w_allp{1},1);
alpha=0.25;
idx=alpha*nperm;

tag_info=cell(nbses,1);
for sess=1:nbses
    
    w=w_alltr{sess};
    wp=w_allp{sess};
    N=length(w);
    %%
    tag=zeros(N,1);
    for n=1:N
        
        sorted=sort(wp(:,n));
        lb=sorted(idx);         % lower bound
        ub=sorted(end-idx+1);   % upper bound
        tag(n)=or(w(n)<lb,w(n)>ub);
        
    end 
    tag_info{sess}=int8(tag);
    
end
%%

perc_info=sum(cell2mat(tag_info))./numel(cell2mat(tag_info));

%%
if saveres==1
    
    savename=['tag_info_s+c_', sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/weight/tags/';    
    save([savefile,savename],'tag_info')
    
end

