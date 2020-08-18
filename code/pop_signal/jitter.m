% reconstruct with random weights an spike timing (jitter in time)

close all
clear all 
clc

place=1;

if place==0
    showfig=0;
    saveres=1;
else
    saveres=0;
    showfig=1;
    factor=100;
end

nperm=1000;

ba=1;
period=2;
T=400;
 
%%
d_range=[5,10,20,40,80,200];
d=d_range(6);
timevec=1:d:T;
nd=length(timevec);

namebeh={'same','different'};
namea={'V1','V4'};
namep={'tar','test'};

task=['reconstruct with jitter of ',sprintf('%1.0i',d),' ms'];
disp(task)

L=100;                              % length of the kernel                                  % exponential kernel for convolution
tau_vec=0:L;                        % support
lambda=1/20;                        % time constant
kernel=exp(-lambda.*tau_vec); 

%% load weights and spikes

addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/') 
end

loadname=['weight_bac_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',T)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_', sprintf('%1.0i',T)];
load(loadname2,'spikes_inm','spikes_cnm')

ncv=size(spikes_cnm,2);
nbses=length(w_all);

%% reconstruct the signal with permutation

tic

T=size(spikes_cnm{1},3);

xp_diff=zeros(nbses,nperm,T);
xp_same=zeros(nbses,nperm,T);
    
parfor sess=1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spikes=cat(2,spike_diff,spike_same); 
    
    weights=squeeze(w_all{sess});
    
    x_perms=zeros(nperm,2,T);
    for perm=1:nperm
           
        for c=1:2
            rt_vec=[];
            for t=1:nd
                
                idx1=timevec(t);                                              
                
                rp=randperm(d);
                rt=idx1-1+rp;                                               % permute spike timing within the interval of d ms
                rt_vec=cat(2,rt_vec,rt);                                    % collect across intervals
            end
            spikes(:,c)=cellfun(@(x) x(:,:,rt_vec),spikes(:,c),'UniformOutput',false); % redefine spikes
        end
       %%
        %%%%%%%%%%%%%%%%%%%%%% compute LDR %%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [x_cv] = reconstruct_fun(weights,spikes,kernel);
        x_perms(perm,:,:) = x_cv;                    
         
    end
    
    xp_diff(sess,:,:)=x_perms(:,1,:);
    xp_same(sess,:,:)=x_perms(:,2,:);
    
end
toc

%% plot

if showfig==1
    
    y1=squeeze(nanmean(xp_same))*factor;                                        % mean across sessions
    y2=squeeze(nanmean(xp_diff))*factor;
    
    figure()
    hold on
    for perm=1:nperm
        plot(y1(perm,:),'g')
        plot(y2(perm,:),'b')
    end
    hold off
    box off
    
    %ylim([-0.5,0.5])
    
    xlabel('time (ms)')
    ylabel('signal from permuted')
    
end
%}


%% save

if saveres==1
        
    savename=['jitter_window_',sprintf('%1.0i',d)];
    savefile='/home/veronika/synced/transfer_result/signal/permute_timing/';
    save([savefile,savename],'xp_diff','xp_same','d_range')
    %clear all 
    
end
    
%}


