% compute weights for two classification problems, S+C and C; with sampling

clear all
close all
clc
  
place=1;
saveres=0;
showfig=1;

ba=1;
K=400;

%%                                                                          
nfold=10;                                                                           % number of folds for the cross-validation (search of the best regularization param.)
Cvec=[0.0012,0.00135,0.0015,0.002,0.005,0.01,0.05,0.1,0.5];

namei={'S+C','C'};
beh={[1,3],[1,2]};

addpath('/home/veronika/synced/transfer_result/input/bootstrapped/');
if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/');
else
    addpath('/home/veronika/transfer/code/function/');
end

loadname=['sc_bigboot_', sprintf('%1.0i',K)];
load(loadname)

%% compute weights of the linear SVM with monte-carlo cv

task=['compute weight comparison K= ' , sprintf('%1.0i',K)];
disp(task)

nbses=size(sc_bootstrap,2);
ncv=size(sc_bootstrap,3);

tic

w_comp=cell(2,nbses);
for type=1:2
    
    sc_diff=squeeze(sc_bootstrap(beh{type}(1),:,:));
    sc_same=squeeze(sc_bootstrap(beh{type}(2),:,:));
    
    parfor sess=1:nbses                                                                  % [sessions]
        warning('off','all');
         
        N=size(sc_same{sess,1},2);
        weight_cv=zeros(ncv,N);
        
        for cv= 1:ncv                                                               % [cross-validation for splits into training and testing]
            
            s1=sc_same{sess,cv}; % decision same, training
            s2=sc_diff{sess,cv}; % -1 decision different training
            
            %%%%%%%%%%%%% get weights, bias and balanced accuracy
            
            [wtilde] = weights2_fun(s1,s2,nfold,Cvec);
            weight_cv(cv,:)=wtilde;
            
        end
        
        w_comp{type,sess}=weight_cv;
        
    end
end

toc


%% save result

if saveres==1
    
    savename=['weight_comp2_',sprintf('%1.0i',K)];
    savefile='/home/veronika/transfer/result/weight/weight_comp2/';    
    save([savefile,savename],'w_comp')
    
end
