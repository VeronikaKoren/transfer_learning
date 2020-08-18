% precompute spike counts and spike trains for training and testing

clear all
close all
clc
                                                                  
saveres=0;

ba=1;     
K=500;                                                              % length of the time window
                                                 
beh=[1,2,3];
nboot=100;                                                          % nb of bootstraps

%%
                                                                            
start=500;                                                           % start time
                                                                     
namea={'V1','V4'};
ending={'_all','_lay'};
namevar={'spikes_test';'spikes_testV4_lay'};

dname=['/home/veronika/v1v4/data/',namea{ba},ending{ba},'/'];
    
addpath(dname)
fname=dir([dname filesep '*.mat']);
nbses=length(fname);

task=['get spike counts bootstrapped K = ',sprintf('%1.0i',K)] ;
disp(task)

%% get boostrapped samples

sc_bootstrap=cell(3,nbses,nboot);

for sess=1:nbses                                                                   % loop across sessions
    
    idxbeh=1:3;
    s=load([dname filesep fname(sess).name],namevar{ba});                                  % load spike trains
    %disp(sess)
    
    x=s.(namevar{ba});
    x_col=cellfun(@(x,y,z) cat(2,x,y,z),x(beh,1),x(beh,2),x(beh,3),'UniformOutput',false); % concatenate layers
    x_time=cellfun(@(x) x(:,:,start:start+K-1),x_col,'UniformOutput',false);               % take the time window
    N=size(x_time{1},2);
    
    sc_all=cellfun(@(x) sum(x,3), x_time, 'UniformOutput', false);
    
    Jall=cellfun(@(x) size(x,1), sc_all);
    [val,idxmax]=max(Jall);
    
    Jmax=Jall(idxmax);
    idxbeh(idxmax)=[];
    
    sc_bootstrap(idxmax,sess,:)=repmat(sc_all(idxmax),nboot,1);                           % condition with least trials, replicate nboot-times
    Jall(idxmax)=[];
    %% bootstrap
    
    %{                                                                                       % bootstrap the other two conditions: 
    for c=1:2
        
        x=sc_all{idxbeh(c)};
        
        for b=1:nboot
            rp=randi(Jall(c),Jmax,1);                                                    % random permutation with replacement, Jmax trials
            xboot=x(rp,:);                                                        
            
            
            sc_bootstrap{idxbeh(c),sess,b}=xboot;
        end
    end
    %}
    
end
    
%% save result
    
if saveres==1
    savename=['sc_bigboot_', sprintf('%1.0i',K)];
    savefile='/home/veronika/transfer/result/input/bootstrapped/';
    save([savefile,savename],'sc_bootstrap')
end
    
    

