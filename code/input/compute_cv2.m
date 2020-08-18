% compute CV2 in single neurons in conditions CNM and INM and test every neuron

clear all
close all
clc
                                                                  
saveres=0;

ba=1; 
period=1;
cond=[1,2];                       % conditions CNM and INM                                                 
    
%%

if period==1
    time_vec=200:599;                                                                            
else
    time_vec=500:899;
end
T=length(time_vec);
                                          
namep={'target','test'};
namea={'V1','V4'};
cond_all={'CNM','INM','CM'};
namec=cond_all(cond);
ending={'_all','_lay'};

variables={'spikes_tar','spikes_test';'spikes_tarV4_lay','spikes_testV4_lay'};

%% compute mean pm SEM of CV2

task=['compute cv2 in ', namea{ba}, ' during ',namep{period}];
display(task)

dname=['/home/veronika/v1v4/data/',namea{ba},ending{ba},'/'];
addpath(dname)
fname=dir([dname filesep '*.mat']);

addpath('/home/veronika/Dropbox/matlabfun/nancorr/')
%%
nbses=length(fname);
cvar=variables{ba,period};

dt=1/1000;

cv_all=[];

for sess=1:nbses                                                               % loop across sessions
    
    s=load([dname filesep fname(sess).name],cvar);                                      % load spike trains
    
    x=s.(cvar)(cond,:);                                                                % take 2 conditions
    x_col=cellfun(@(x,y,z) cat(2,x,y,z),x(:,1),x(:,2),x(:,3),'UniformOutput',false);    % concatenate layers
    
    %% get spike counts for training and spike trains for testing with monte-carlo cv
    
    N=size(x_col{1},2);
    J=cellfun(@(x) size(x,1),x_col);                                                   			% number of trials
    
    x_time=cellfun(@(x) x(:,:,time_vec),x_col,'UniformOutput',false);
    x_perm=cellfun(@(x) x(randperm(size(x,1)),:,:),x_time, 'UniformOutput',false);                            % permute the order of trials
    x_cut=cellfun(@(x) x(1:min(J),:,:), x_perm, 'UniformOutput',false);                                   % take same number of trials
    
    %%
    cv_n=cell(N,2);
    
    for c=1:2   
        for n=1:N
            
            cv_m=[];
            
            for j=1:J(c)
                
                xs=squeeze(x_time{c}(j,n,:));
                st=find(xs);
                isi=st(2:end)-st(1:end-1);
                cv2=2*abs(isi(2:end)-isi(1:end-1))./(isi(2:end)+isi(1:end-1));
                
                if isempty(cv2)==0
                    cv_m=cat(1,cv_m,cv2);           % concatenate cv2 across trials
                end
            end
            
            cv_n{n,c}=cv_m;            
            
        end
    end
    
    cv_all=cat(1,cv_all,cv_n);
      
end

%% mean and SEM

cv_mean=cellfun(@(x) nanmean(x), cv_all);
cv_sem=cellfun(@(x) nanstd(x)./sqrt(size(x,1)), cv_all);


%% remove NaN

nidxs=[];
for c=1:2
    
    nidx=find(isnan(cv_mean(:,c)));
    nidxs=cat(1,nidxs,nidx);
    
end

cv_mean(nidxs,:)=[];
cv_sem(nidxs,:)=[];

R=corr(cv_mean(:,1),cv_mean(:,2));

%% sort

[val,idx]=sort(cv_mean(:,1));
order=flip(idx);

cv_sorted=cv_mean(order,:);
sem_sorted=cv_sem(order,:);

%% save result

if saveres==1
    savename=['cv2_',namea{ba},'_', namep{period}];
    savefile='/home/veronika/synced/transfer_result/basic_stat/cv2/';
    save([savefile,savename],'cv_sorted','sem_sorted','R','namec')
end

%%

