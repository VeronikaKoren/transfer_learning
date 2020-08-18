% compute the performance of the linear SVM on sign-specific model: extract
% the permformance of a specific group by replacing the input of the other
% group with activity from a randomly selected trial

close all
clear all
clc
format long

%%%%%%%%%%%

place=0;                                                                        % 0 for the server, 1 for the office computer

ba=1;                                                                           % brain area: 1 for V1, 2 for V4
period=2;

saveres=1;                                                                      % save result?

nperm=100;                                                                        % number of permutations 
ncv=100;                                                                         % number of cross-validations for splits into training and validation set 

%%
info_case=1;

Cvec=[0.0012,0.00135,0.0015,0.002,0.005,0.01, 0.05, 0.1,0.5];                       % range of tested regularization parameters
ratio_train_val=0.8;                                                            % ratio of training/validation data
nfold=5;                                                                       % number of folds for computing the regularization param

start_vec=[200,500];
start=start_vec(period);                                                      % start of the time window
K=500; 

%%
namea={'V1','V4'};
namep={'target','test'};
names={'minus','plus'};
namei={'s+c','c'};

addpath('/home/veronika/synced/transfer_result/input/spike_train/')
addpath('/home/veronika/synced/transfer_result/weight/weight_alltr/')
if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else    
    addpath('/home/veronika/transfer/code/function/')
end

%% load

% load spike counts
loadname=['spike_train_s+c_',namea{ba},'_',namep{period}];
load(loadname)
x_time=cellfun(@(x) x(:,:,start:start+K-1),spiketrain, 'UniformOutput', false);
sc_all=cellfun(@(x) sum(x,3), x_time,'UniformOutput', false);
nbses=size(sc_all,1);

% load weights
loadname2=['weight_alltr_',namei{info_case},'_',namea{ba} ,'_', sprintf('%1.0i',K)];
load(loadname2,'w_alltr')

signw=cellfun(@sign, w_alltr, 'UniformOutput', false);

tag=[-1,1];
%% classify

display(['computing svm sign remove info ', namea{ba},' ',namep{period}])

tic
bac_sign=cell(nbses,1);

parfor sess=1:nbses
   
    sc_sess=sc_all(sess,:);                 % 2 conditions
    sc_one=cat(1,sc_sess{1},sc_sess{2});    % concatenated conditions
    
    N=size(sc_one,2);
    J_tot=size(sc_one,1);
    J=cellfun(@(x) size(x,1),sc_sess);
    
    bac_s=zeros(2,nperm);
    for sgn=1:2

        idx_keep=find(signw{sess}==tag(sgn));
        
        for perm=1:nperm
           
            sc_perm=sc_one(randperm(J_tot),:,:);                                    % random permutation of the trial order
            sc_perm(:,idx_keep)=sc_one(:,idx_keep);                                 % selected group has correct order
            s1=sc_perm(1:J(1),:);
            s2=sc_perm(J(1)+1:end,:);
            [bac] = svm_simple_fun(s1,s2,ratio_train_val,ncv,nfold,Cvec);
            
            bac_s(sgn,perm)=bac;                                                    % balanced accuracy
        end                                                                         
         
    end
    
    bac_sign{sess}=bac_s;
end

toc

%%

bac_minus=cellfun(@(x) x(1,:),bac_sign,'UniformOutput', false);
bac_plus=cellfun(@(x) x(2,:),bac_sign,'UniformOutput', false);

%% save results

if saveres==1
    address='/home/veronika/synced/transfer_result/weight/bac/';
    filename=['bac_sign_', namea{ba},namep{period},'_',sprintf('%1.0i',K)];
    save([address, filename],'bac_minus','bac_plus','K','start')
end

