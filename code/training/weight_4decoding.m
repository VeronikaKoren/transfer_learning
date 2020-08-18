% compute feature weights of the linear SVM in sessions
% use pre-computed data for monte-carlo cross-validation 
% REGULAR MODEL

clear all
close all
clc
  
place=1;
saveres=0;
showfig=1;

ba=1;
period=2;
T=400;

info_case=1;    % 1 for stim+choice and 2 for choice
%%                                                                          

namea={'V1','V4'};
namep={'tar','test'};
namec={'','c'};

nfold=10;                                                                           % number of folds for the cross-validation (search of the best regularization param.)
Cvec=[0.0012,0.00135,0.0015,0.002,0.005,0.01,0.05,0.1,0.5];

addpath('/home/veronika/synced/transfer_result/input/');
if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/');
else
    addpath('/home/veronika/transfer/code/function/');
end

loadname=['input_',namec{info_case},namea{ba},'_',namep{period},'_', sprintf('%1.0i',T)];
load(loadname)

nbses=size(count_cnm,1);
ncv=size(count_cnm,2);

if info_case==1
    % training data : stim+choice
    train_diff=count_cnm;
    train_same=count_cm;
    
    % testing data: choice
    test_diff=cellfun(@(x) sum(x,3),spikes_cnm,'UniformOutput', false);
    test_same=cellfun(@(x) sum(x,3),spikes_inm,'UniformOutput', false);
    
else
    
    % training data: choice
    train_diff=count_cnm;
    train_same=count_inm;
    
    % testing data: choice
    test_diff=cellfun(@(x) sum(x,3),spikes_cnm,'UniformOutput', false);
    test_same=cellfun(@(x) sum(x,3),spikes_inm,'UniformOutput', false);
    
end

%% compute weights of the linear SVM with monte-carlo cv

task=['compute weights in ', namea{ba}, namep{period}, ' T=', sprintf('%1.0i',T)];
disp(task)

tic

w_all=cell(nbses,1);

parfor sess=1:nbses                                                                  % [sessions]
    warning('off','all');
    
    s1=train_same{sess}; % decision same, training
    t1=test_same{sess}; % decision same, testing
    
    N=size(s1,2);
    weight_cv=zeros(ncv,N);
    bac_cv=zeros(ncv,1);
    
    for cv= 1:ncv                                                                 % [cross-validation for splits into training and testing]
        
        s2=train_diff{sess,cv}; % -1 decision different training 
        t2=test_diff{sess,cv}; % -1 decision different, testing
        
        %%%%%%%%%%%%% get weights, bias and balanced accuracy
        
        [wtilde] = weights2_fun(s1,s2,nfold,Cvec,t1,t2); 
        
        weight_cv(cv,:)=wtilde;
       
    end
    
    w_all{sess}=weight_cv;
                                                                % collect across sessions
end

toc


%% show plot

if showfig==1
    figure('Units','centimeters', 'Position', [0,0,24,24])
    for sess=1:nbses
        subplot(4,5,sess)
        boxplot(w_all{sess})
        line([0,size(w_all{sess},2)+1],[0,0],'color','k','linestyle','--')
        box off
    end
    axes


    h1 = xlabel ('neuron index','units','normalized','Position',[0.5,-0.08,0]);
    h2 = ylabel ('decoding weight','units','normalized','Position',[-0.12,0.5,0]);
    h3=title('decoding weights in sessions, distr. over cv');
    
    set(gca,'Visible','off')
    set(h1,'visible','on')
    set(h2,'visible','on')
    set(h3,'visible','on')

    
end
%% save result

if saveres==1
    
    savename=['weight_choice_',namea{ba},'_',namep{period}, '_',sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/weight/weight_decoding/';    
    save([savefile,savename],'w_all')
    
end


