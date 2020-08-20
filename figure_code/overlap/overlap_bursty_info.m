
clear all
close all
clc
                                                                  
K=400;

%%

namea={'V1','V4'};
namet={'bursty','informative'};

%%
addpath('/home/veronika/transfer_learning/result/weight/tags/')
addpath('/home/veronika/transfer_learning/result/weight/weight_alltr/')

loadname=['tag_info_s+c_',sprintf('%1.0i',K)];
load(loadname)

loadname2=['tag_bursty_',sprintf('%1.0i',K)];
load(loadname2)

%%
info=cell2mat(tag_info);
burst=cell2mat(bursty);

idx_burst=find(burst);
idx_info=find(info);

%%
both=cell(2,1);

both{1}=idx_burst;
both{2}=idx_info;

[Nmore,idx_more]=max([numel(idx_burst),numel(idx_info)]);
idx_less=1:2;
idx_less(idx_more)=[];

%%
overlap=0;
for i=1:Nmore
   
    if find(both{idx_more}(i)==both{idx_less},1)
        overlap=overlap+1;
    end
end

perc1=overlap/numel(both{idx_less})*100;
display([sprintf('%0.4f',perc1),' percent of ', namet{idx_less},' neurons is ',namet{idx_more}])

perc2=overlap/numel(both{idx_more})*100;
display([sprintf('%0.4f',perc2),' percent of ', namet{idx_more},' neurons is ',namet{idx_less}])

