% compute feature weights of the linear SVM in sessions using all trials 
% to determine informative and uninformative neurons
% permuted class labels

clear all
close all
clc
  
place=1; 

if place==1
    saveres=0;
    showfig=1;
else
    saveres=1;
    showfig=0;
end


nperm=1000;

cond=2;
ba=1;
period=2;
T=400;

%%                                                                          
namea={'V1','V4'};
namec={'s+c','c'};

nfold=10;                                                                           % number of folds for the cross-validation (search of the best regularization param.)
Cvec=[0.0012,0.00135,0.0015,0.002,0.005,0.01,0.05,0.1,0.5];

addpath('/home/veronika/synced/transfer_result/input/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/');
else
    addpath('/home/veronika/transfer/code/function/');
end

loadname=['sc_',namec{cond},'_',namea{ba},'_',sprintf('%1.0i',T)];
load(loadname)
nbses=size(sc_all,1);

%% compute weights of the linear SVM with monte-carlo cv

task=['compute weights alltr ',namec{cond},'_', namea{ba},' T=', sprintf('%1.0i',T)];
disp(task)

tic
w_allp=cell(nbses,1);
parfor sess=1:nbses                                                                  % [sessions]
    warning('off','all');
    
    s1=sc_all{sess,1}; % condition cm
    s2=sc_all{sess,2}; % condition cnm
    N=size(s1,2);
    J=[size(s1,1),size(s2,1)];
    
    
    wperm=zeros(nperm,N);
    y_correct=cat(1,ones(J(1),1).*(-1),ones(J(2),1));
    
    for perm=1:nperm
        y_train=y_correct(randperm(sum(J)));
        [wtilde] = weights_label_fun(s1,s2,nfold,Cvec,y_train);
        wperm(perm,:)=wtilde;    
    end
    
    w_allp{sess}=wperm;             % collect across sessions
    
                                                                                
end
toc
%%
if showfig==1
    
    figure()
    hold on
    for perm=1:10
        plot(w_allp{1}(perm,:),'k')
    end
    
end
%% save result

if saveres==1
    
    savename=['weight_allp_',namec{cond},'_',namea{ba},'_', sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/weight/weight_alltr/';    
    save([savefile,savename],'w_allp')
    %clear all
end
