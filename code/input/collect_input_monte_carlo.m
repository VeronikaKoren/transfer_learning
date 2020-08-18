% precompute spike counts and spike trains for training and reconstruction

clear all
close all
clc
                                                                  
saveres=0;
ba=1;     
period=2;

T=400;                                                                  % length of the time window
%%
start_vec=[200,500];                                                                              
start=start_vec(period);                                                 % start time
                                                                       
namep={'tar','test'};
namea={'V1','V4'};
ending={'_all','_lay'};
namevar={'spikes_tar','spikes_tarV4_lay';'spikes_test','spikes_testV4_lay'};

ratio=0.5;                                                                   % ration nb trials for train and for test
ncv=100;                                                                     % nb cross-validations

%% 
 
dname=['/home/veronika/v1v4/data/',namea{ba},ending{ba},'/'];
    
addpath(dname)
fname=dir([dname filesep '*.mat']);
nbses=length(fname);

task=['get input with monte-carlo cv for condition CNM in ',namea{ba}, ' ' namep{period}];
disp(task)

%% get data for training and testing

count_cnm=cell(nbses,ncv);
spikes_cnm=cell(nbses,ncv);
count_cm=cell(nbses,1);
spikes_inm=cell(nbses,1);

for sess=1:length(fname)                                                               % loop across sessions
    
    s=load([dname filesep fname(sess).name],namevar{period,ba});  
    %s=load([dname filesep fname(sess).name]);
    %disp(sess)
    
    x=s.(namevar{period,ba});
    x_col=cellfun(@(x,y,z) cat(2,x,y,z),x(1:3,1),x(1:3,2),x(1:3,3),'UniformOutput',false);    % concatenate layers
    x_time=cellfun(@(x) x(:,:,start:start+T-1),x_col,'UniformOutput',false);            % take the time window
    
    
    %% get spike counts for training and spike trains for testing in condition CNM with monte-carlo cv
    Jcnm=size(x_time{1},1);
    idx_cut=round(Jcnm/2);                                                              % index for dividing the data in CNM into training and testing
    
    for cv=1:ncv
        random_order=randperm(Jcnm);
        cnm=x_time{1}(random_order,:,:);                                                % permute the order of trials in condition CNM
        count_cnm{sess,cv}=single(squeeze(sum(cnm(1:idx_cut,:,:),3)));                  % counts CNM (for training)
        spikes_cnm{sess,cv}=int8(cnm(idx_cut+1:end,:,:));                             % spikes CNM (for testing)
    end
    
    %% get counts for training (condition CM) and spikes for testing (condition INM)
    
    count_cm{sess}=single(squeeze(squeeze(sum(x_time{3},3))));                          % counts condition CM
    spikes_inm{sess}=int8(x_time{2});                                                 % spikes coundition INM
    
end
    
%% save result
    
if saveres==1
    savename=['input_',namea{ba},'_',namep{period},'_', sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/input/transfer_mc/';
    save([savefile,savename],'count_cnm','count_cm','spikes_cnm','spikes_inm')
end
    
    

