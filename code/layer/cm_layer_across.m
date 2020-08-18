% cross-mapping across layers 

close all
clear all
clc

place=1;
saveres=0;                                                            % save result?

condition=3;
K=350;

%%
beh_vec=[1,2,3];
beh=beh_vec(condition);

ba=1;
period=2; 
                                                                      % number of time steps
start_vec=[200,500]+140;                                              % beginning of the time window for the target (200) and the test stimulus (500) 
start=start_vec(period);

D=50;                                                                 % number of dimensions (a.k. max delay)

support=-10:10;
sigma_u=length(support)/5;
u=exp(-(support).^2/(2*sigma_u));                                     % gaussian kernel
u=u./sum(u);    

namea={'V1','V4'};
namep={'target','test'};
namec={'cnm','inm','cm'};   
nameg={'SG','G','IG'};

display(['cross-mapping across layers ',namec{condition}])

%% add path

addpath('/home/veronika/synced/transfer_result/input/spike_train/');
addpath('/home/veronika/synced/transfer_result/weight/tags/');
addpath('/home/veronika/synced/transfer_result/signal/layer/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load

loadname=['spike_train_all_',namea{ba},'_',namep{period}];   % all spike trains
load(loadname);

ssingle=cellfun(@(x) single(x),spiketrain(:,condition), 'un', 0);                                         % take 1 condition
stime=cellfun(@(x) x(:,:,start:start+K-1), ssingle,'UniformOutput', false);                          % select time window

loadname3=['ncell_lay_',namea{ba}];
load(loadname3)
nbses=size(stime,1);
%%                                                                          
idx1=[1,1,2];
idx2=[2,3,3];

labels=cell(3,2);

for c=1:3
    str1=nameg{idx2(c)};
    str2=nameg{idx1(c)};
    
    labels{c,1}=[str1,' to ', str2];
    labels{c,2}=[str2,' to ', str1];
end

%% cross-mapping

a1=cell(nbses,1);
a2=cell(nbses,1);
a3=cell(nbses,1);

tic

parfor sess = 1:nbses
    
    spikes=cell(3,1);
    for g=1:3
        
        ncl=ncell_lay(sess,:);
        s=cumsum([0,ncl]);
        tagged = s(g) + 1 : s(g+1);
        spikes{g}=stime{sess}(:,tagged,:);
        
    end
    %%  
    
    for c=1:3
        x1=spikes{idx1(c)};
        x2=spikes{idx2(c)};
        [r_cm] = compute_cma_fun(x1,x2,D,u);      % compute cross-mapping from one layer to another
        if c==1
            a1{sess}=r_cm';
        elseif c==2
            a2{sess}=r_cm';
        elseif c==3
            a3{sess}=r_cm';
        end
    end
    
end
    
toc

%% save results

if saveres==1
    address='/home/veronika/synced/transfer_result/pairwise/cm/across/ ';
    filename=['cm_layer_across_',namec{condition},'_',sprintf('%1.0i',K)];
    save([address, filename], 'a1','a2','a3', 'D','u','start','labels')
end

%clear all
%exit
