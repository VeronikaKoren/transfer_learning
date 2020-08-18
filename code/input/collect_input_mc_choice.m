% precompute spike counts and spike trains for training and testing/ choice
% only

clear all
close all
clc
                                                                  
saveres=0;
ba=1;     
period=2;

T=400;                                                                       % length of the time window
%%
start_vec=[200,500];                                                                              
start=start_vec(period);                                                                   % start time

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

task=['get input with monte-carlo cv for conditions CNM/INM in ',namea{ba}, ' ' namep{period}];
disp(task)

%% get data for training and testing

count=cell(nbses,ncv,2);
spikes=cell(nbses,ncv,2);

for sess=1:length(fname)                                                               % loop across sessions
    
    s=load([dname filesep fname(sess).name],namevar{period,ba});  
    
    
    x=s.(namevar{period,ba});
    x_col=cellfun(@(x,y,z) cat(2,x,y,z),x(1:2,1),x(1:2,2),x(1:2,3),'UniformOutput',false);    % concatenate layers
    x_time=cellfun(@(x) x(:,:,start:start+T-1),x_col,'UniformOutput',false);            % take the time window
    
    
    %% get spike counts for training and spike trains for testing (monte-carlo cv)
    J=cellfun(@(x) size(x,1),x_time);
    idx_cut=round(J/2);                                                              % index for dividing the data in CNM into training and testing
    %%
    
    for c=1:2
        for cv=1:ncv
            
            random_order=randperm(J(c));
            xtrain=x_time{c}(random_order,:,:);                                            % permute the order of trials 
            count{sess,cv,c}=single(squeeze(sum(xtrain(1:idx_cut(c),:,:),3)));                % counts (for training)
            spikes{sess,cv,c}=int8(xtrain(idx_cut(c)+1:end,:,:));                             % spikes(for testing)
        end
    end
    
end

%%

count_cnm=squeeze(count(:,:,1));
count_inm=squeeze(count(:,:,2));
spikes_cnm=squeeze(spikes(:,:,1));
spikes_inm=squeeze(spikes(:,:,2));

%% save result
    
if saveres==1
    savename=['input_c',namea{ba},'_',namep{period},'_', sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/input/transfer_mc/';
    save([savefile,savename],'count_cnm','count_inm','spikes_cnm','spikes_inm')
end
    
    

