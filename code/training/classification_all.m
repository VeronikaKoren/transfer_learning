% compute feature weights of the linear SVM in sessions
% use pre-computed data for monte-carlo cross-validation 
% REGULAR MODEL

clear all
close all
clc
  
place=1;
saveres=0;

ba=1;
period=2;
T=400;

info_case=1;    % 1 for stim+choice and 2 for choice
nperm=1000;
%%                                                                          

namea={'V1','V4'};
namep={'tar','test'};
namei={'s+c','c'};

ratio_train_val=0.8;                                                                % ratio of training/validation data
nfold=10;                                                                           % number of folds for the cross-validation (search of the best regularization param.)
Cvec=[0.0012,0.00135,0.0015,0.002,0.005,0.01,0.05,0.1,0.5];
ncv=50;

%% load

addpath('/home/veronika/synced/transfer_result/input/spike_count/');
if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/');
else
    addpath('/home/veronika/transfer/code/function/');
end

loadname=['sc_',namei{info_case},'_',namea{ba},'_',namep{period},'_', sprintf('%1.0i',T)];
load(loadname)

nbses=size(sc_all,1);
%% compute weights of the linear SVM with monte-carlo cv

task=['compute BAC ',namei{info_case},' ', namea{ba}, namep{period}, ' T=', sprintf('%1.0i',T)];
disp(task)

tic

bac_all=zeros(nbses,1);
bac_allp=zeros(nbses,nperm);

parfor sess=1:nbses                                               % sessions
    warning('off','all');
    
    s1=sc_all{sess,1}; % condition CNM
    s2=sc_all{sess,2}; % condition CM
    
    [bac,bacp] = svm_mc_fun(s1,s2,ratio_train_val,ncv,nfold,nperm,Cvec);
    
    bac_all(sess)=bac;
    bac_allp(sess,:)=bacp;                                                           % collect across sessions
end

toc
%%
display(nanmean(bac_all),'bac')
display(nanmean(nanmean(bac_allp)),'bac permuted')

%% save result

if saveres==1
    
    savename=['bac_',namei{info_case},'_',namea{ba},'_',namep{period}, '_',sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/weight/bac/';    
    save([savefile,savename],'bac_all','bac_allp')
    
end


