% cross-mapping within layers 

close all
clear all
clc

place=0;
saveres=0;                                                                  % save result?

info_case=4;
K=500;

%%
beh_info={[1,2,3],1,2,3};
beh=beh_info{info_case};

ba=1;
period=2; 
                                                                      % number of time steps
start_vec=[200,500];                                                  % beginning of the time window for the target (200) and the test stimulus (500) 
start=start_vec(period);

D=50;                                                                 % number of dimensions (a.k. max delay)

support=-10:10;
sigma_u=length(support)/5;
u=exp(-(support).^2/(2*sigma_u));                                     % gaussian kernel
u=u./sum(u);    

namea={'V1','V4'};
namep={'target','test'};
namei={'all','cnm','inm','cm'};   

display(['cross-mapping within layers ',namei{info_case}])

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

if info_case==1
    ssingle=cellfun(@(x,y,z) single(cat(1,x,y,z)),spiketrain(:,1),spiketrain(:,2),spiketrain(:,3), 'un', 0);    % concatenate 3 conditions
else
    ssingle=cellfun(@(x) single(x),spiketrain(:,info_case-1), 'un', 0);                                         % take 1 condition
end

stime=cellfun(@(x) x(:,:,start:start+K-1), ssingle,'UniformOutput', false);                                     % select time window
nbses=size(stime,1);

loadname3=['ncell_lay_',namea{ba}];
load(loadname3)

%%                                                                          

cm_lay=cell(3,1);

tic
for g=1:3
    
    r_sess=cell(nbses,1);
    
    parfor sess = 1:nbses
        
        ncl=ncell_lay(sess,:);
        s=cumsum([0,ncl]);
        tagged = s(g) + 1 : s(g+1);
        
        if length(tagged)>1
            
            spike_train=stime{sess}(:,tagged,:);
            [r_cm] = compute_cm2_fun(spike_train,D,u);           
            r_sess{sess}=r_cm';                                    % collect across sessions
            
        end
    end
    
    cm_lay{g}=cell2mat(r_sess);
        
end

toc

if showfig==1
   figure()
   for g=3
       subplot(3,1,g)
       plot(cm_lay{g}(:,1))
       hold on
       plot(cm_lay{g}(:,2))
       ylim([0,0.7])
   end
end

%% save results

if saveres==1
    address='/home/veronika/synced/transfer_result/pairwise/cm/within/';
    filename=['cm_layer_within_',namei{info_case},'_',sprintf('%1.0i',K)];
    save([address, filename], 'cm_lay', 'D','u','start')
end

%clear all
%exit



