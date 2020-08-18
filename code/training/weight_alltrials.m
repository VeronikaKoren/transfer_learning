% compute feature weights of the linear SVM in sessions using all trials 


clear all
close all
clc
  
saveres=0;
showfig=1;

ba=1;
period=2;
info_case=2;

T=400;
start=500;
%%                                                                          
namea={'V1','V4'};
namep={'target','test'};
namec={'s+c','c'};

task=['compute weights all trials regular ' namec{info_case}];
disp(task)

nfold=10;                                                                           % number of folds for the cross-validation (search of the best regularization param.)
Cvec=[0.0012,0.00135,0.0015,0.002,0.005,0.01,0.05,0.1,0.5];

addpath('/home/veronika/synced/transfer_result/input/spike_count/');
addpath('/home/veronika/Dropbox/transfer/code/function/');

loadname=['sc_',namec{info_case},'_',namea{ba},'_',namep{period},'_',sprintf('%1.0i',T)];
load(loadname)
nbses=size(sc_all,1);

%% compute weights of the linear SVM using all trials

w_alltr=cell(nbses,1);

for sess=1:nbses                                                         
    warning('off','all');
    
    s1=sc_all{sess,1}; % condition 1
    s2=sc_all{sess,2}; % condition 2
    
    [wtilde] = weights2_fun(s1,s2,nfold,Cvec);
    w_alltr{sess}=wtilde;
                                                                  % collect across sessions
end

%% show plot

if showfig==1
    
    figure()
    for sess=1:nbses
        subplot(4,5,sess)
        ksdensity(w_alltr{sess})
    end
end
%% save result

if saveres==1
    
    savename=['weight_alltr_',namec{info_case},'_',namea{ba},'_', sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/weight/weight_alltr/';    
    save([savefile,savename],'w_alltr')
    
end
