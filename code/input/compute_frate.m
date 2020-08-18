
clear all
close all
clc
                                                                  
saveres=0;
showfig=0;

ba=1; 
period=2;
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

%% compute the mean pm SEM firing rate of single neurons

task=['compute frate in ', namea{ba}, ' during ',namep{period}];
display(task)

dname=['/home/veronika/v1v4/data/',namea{ba},ending{ba},'/'];
addpath(dname)
fname=dir([dname filesep '*.mat']);

nbses=length(fname);
dt=1/1000;
cvar=variables{ba,period};

frate=cell(nbses,1);
sem=cell(nbses,1);
J_all=zeros(nbses,2);

for sess=1:nbses                                                            
    
    s=load([dname filesep fname(sess).name],cvar);                                      % load spike trains
    
    x=s.(cvar)(cond,:);                                                                % take 2 conditions
    x_col=cellfun(@(x,y,z) cat(2,x,y,z),x(:,1),x(:,2),x(:,3),'UniformOutput',false);    % concatenate layers
    
    %% get spike counts for training and spike trains for testing with monte-carlo cv
    
    N=size(x_col{1},2);
    J=cellfun(@(x) size(x,1),x_col);                                                   % number of trials
    
    x_time=cellfun(@(x) x(:,:,time_vec),x_col,'UniformOutput',false);
    x_perm=cellfun(@(x) x(randperm(size(x,1)),:,:),x_time, 'UniformOutput',false);                            % permute the order of trials
    x_cut=cellfun(@(x) x(1:min(J),:,:), x_perm, 'UniformOutput',false);                                   % take same number of trials
    
    x_count=cellfun(@(x) squeeze(mean(x,3))./dt, x_cut, 'UniformOutput',false);                                   % spike count
    x_std=cellfun(@(x) std(x)./sqrt(min(J)), x_count, 'UniformOutput',false);
    x_mean=cellfun(@(x) mean(x), x_count, 'UniformOutput',false);
    
    J_all(sess,:)=J;
    
    frate{sess}=double(cell2mat(x_mean)');  % collect across sessions
    sem{sess}=double(cell2mat(x_std)');
    
end

%% put in matrix and sort

fr_all=cell2mat(frate);
sem_all=cell2mat(sem);

[~,idx]=sort(fr_all(:,1));
order=flip(idx);

fr_sorted=fr_all(order,:);
sem_sorted=sem_all(order,:);

%% plot 

if showfig==1
    
    col={'b','g'};
    x=1:length(idx);
    
    figure()
    for i=1:2
        y1=fr_sorted(:,i) - sem_sorted(:,i);
        y2=fr_sorted(:,i) + sem_sorted(:,i);
        
        hold on
        patch([x fliplr(x)], [y1' fliplr(y2')], col{i},'FaceAlpha',0.3,'EdgeColor',col{i})
        
    end

end
%%

display(mean(J_all),'mean nb trial');

%% save result

if saveres==1
    savename=['frate_',namea{ba},'_', namep{period}];
    savefile='/home/veronika/synced/transfer_result/basic_stat/frate/';
    save([savefile,savename],'fr_sorted','sem_sorted','namec')
end

%%

