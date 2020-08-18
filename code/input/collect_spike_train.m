
close all
clear all
clc

saveres=0;                                                                          % save result?

info_case=1;                                                                    % 1 for S+C, 2 for C
beh_info={[1,3],[1,2]};

beh=beh_info{info_case};
%%

ba=1;
period=2;

start_vec=[200,500];                                                             % beginning of the time window for the target (200) and the test stimulus (500) 
start=start_vec(period);

ending={'_all','_lay'};
namea={'V1','V4'};
namep={'target','test'};
namei={'s+c','c','all'};
variables={'spikes_tar','spikes_test';'spikes_tarV4_lay','spikes_testV4_lay'};
                                                                     
display(['collect spike trains ',namei{info_case},' ',namea{ba},' ', namep{period}])

%%

addpath('/home/veronika/v1v4/data/')
dname=['/home/veronika/v1v4/data/',namea{ba},ending{ba},'/'];    
fname=dir([dname filesep '*.mat']);                                             % name of the data file                                                                                                                                                                                            
cvar=variables{ba,period};

%%                                                                         
spiketrain=[];

for sess = 1:length(fname)
    
    display(sess,'session')
    
    s=load([dname filesep fname(sess).name],cvar);                                                                                      % load data in a session
    %N=size(s.(cvar){1,1},2)+size(s.(cvar){1,2},2)+size(s.(cvar){1,3},2);
    
    s_beh=s.(cvar)(beh,:);
    
    s_catlay=cellfun(@(x,y,z) cat(2,x,y,z), s_beh(:,1), s_beh(:,2), s_beh(:,3), 'UniformOutput', false);
    s_int=cellfun(@(x) int8(x), s_catlay', 'UniformOutput', false);
    spiketrain=cat(1,spiketrain,s_int);
    
end

%% save results

if saveres==1
    address='/home/veronika/synced/transfer_result/input/spike_train/';
    filename=['spike_train_',namei{info_case},'_',namea{ba},'_',namep{period}];
    save([address, filename], 'spiketrain')
   
end


