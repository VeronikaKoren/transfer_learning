
clear all
close all
clc
                                                                  
saveres=1;

ba=1;     
period=2;
T=400;
 
%%
namep={'tar','test'};
namea={'V1','V4'};

addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/')

loadname=['weight_alltr_',namea{ba},'_', sprintf('%1.0i',T)];
load(loadname)

%%
plus_alltr=cellfun(@(x) x>0,w_alltr,'UniformOutput', false);

if saveres==1
    savename=['tag_sign_',namea{ba},'_', sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/weight/weight_alltr/';
    save([savefile,savename],'plus_alltr')
end
