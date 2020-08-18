% pre-compute spike counts  

clear all
close all
clc
                                                                  
saveres=0;

ba=1;     
period=2;

info_case=2;

%%
                                                                            
start=500;                                                                   % start time
T=400;                                                                       % length of the time window

behaviors={[1,3],[1,2]};
beh=behaviors{info_case};

namea={'V1','V4'};
namep={'target','test'};
ending={'_all','_lay'};
namevar={'spikes_tar','spikes_tarV4_lay';'spikes_test','spikes_testV4_lay'};
namei={'s+c','c'};

%% 
 
dname=['/home/veronika/v1v4/data/',namea{ba},ending{ba},'/'];    
addpath(dname)
fname=dir([dname filesep '*.mat']);
nbses=length(fname);

task=['get spike counts all trials in ',namei{info_case},' ' namea{ba}];
disp(task)

%% get spike counts in all conditions

sc_all=cell(nbses,length(beh));

for sess=1:nbses                                                         % loop across sessions
    
    s=load([dname filesep fname(sess).name],namevar{period,ba});                % load spike trains
    %disp(sess)
    
    x=s.(namevar{period,ba});
    x_col=cellfun(@(x,y,z) cat(2,x,y,z),x(beh,1),x(beh,2),x(beh,3),'UniformOutput',false); % concatenate layers
    x_time=cellfun(@(x) x(:,:,start:start+T-1),x_col,'UniformOutput',false);               % take the time window
    %N=size(x_time{1},2)
    
    sc_all(sess,:)=cellfun(@(x) sum(x,3), x_time, 'UniformOutput', false);
    
end
    
%% save result

if saveres==1
    savename=['sc_',namei{info_case},'_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
    savefile='/home/veronika/synced/transfer_result/input/spike_count/';
    save([savefile,savename],'sc_all')
end
    
    

