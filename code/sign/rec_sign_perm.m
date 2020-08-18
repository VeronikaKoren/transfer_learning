% reconstruct signal sign

clear all 
close 
clc
format short

place=1;

if place==0
    saveres=1;
    showfig=0;
    
else
    saveres=1;
    showfig=1;
    factor=100;
    
end
 
nperm=2;

ba=1;
period=2;
tau=20;                                                                 
T=500;
%%

namea={'V1','V4'};
namep={'tar','test'};

task='compute sign permuted ';
disp(task)

% exponential kernel for convolution
L=100;                                                                       % length of the kernel                                  
tau_vec=0:L;                                                                 % support
lambda=1/tau;                                                                % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

if place==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end
addpath('/home/veronika/synced/transfer_result/input/')                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/');

%% load results

loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',T)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',T)];
load(loadname2,'spikes_inm','spikes_cnm')
      

nbses=length(w_all);
ncv=size(w_all{1},1);

%%

xp_minus=zeros(nbses,nperm,2,T);
xp_plus=zeros(nbses,nperm,2,T);


tic
parfor sess = 1:nbses
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);             % concatenate conditions
    
    J=[size(spike_diff{sess,1},1), size(spike_same{sess,1},1)];
    N=size(spike_diff{1},2);
    %%
    w_regular=squeeze(w_all{sess});
    w_abs=max(max(abs(w_regular)));
    
    x_min=zeros(nperm,2,T);
    x_pl=zeros(nperm,2,T);
    
    for perm=1:nperm                                                                             % permutation cycles
        
        w_random=-w_abs + (w_abs + w_abs).*rand(1,N);
        
        spike_perm=cellfun(@(x) x(randperm(sum(J)),:,:),spike_one, 'UniformOutput',false);        % permute class label for spike trains
        spikes=cell(ncv,2);
        spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spike_perm, 'UniformOutput',false);
        spikes(:,2)=cellfun(@(x) x(J(1)+1:sum(J),:,:), spike_perm, 'UniformOutput',false);
    
        for r = 1:2
                
            weights=repmat(w_random,ncv,1);
            [x_cv] = reconstruct_fun(weights,spikes,kernel);
            
            if r==1
                x_min(perm,:,:)=x_cv;                                                 % collect results across permutation
            elseif r==2
                x_pl(perm,:,:)=x_cv;
            end
            
        end
    end
    % perm
    
    xp_minus(sess,:,:,:)=x_min;                                                % collect results across sessions
    xp_plus(sess,:,:,:)=x_pl;
    
end

toc


%% plot

if showfig==1
    
    %x_show=squeeze(xp_minus(sess,:,2,:)-xp_minus(sess,:,1,:)).*factor;
    x_show=squeeze(xp_plus(1,:,2,:)-xp_plus(1,:,1,:)).*factor;
    
    y1=min(x_show);
    y2=max(x_show);
    
    figure()
    hold on
    plot(y1,'k')
    plot(y2,'k')
    
    hold off
    box off
    
    ylim([-0.5,0.5])
    
    xlabel('time (ms)')
    ylabel('signal from permuted')
    
end
%}

%%

xp_sign=cell(2,1);   % put in a cell
xp_sign{1}=xp_minus; 
xp_sign{2}=xp_plus;

%% save

if saveres==1
    
    savename=['sign_perm_', sprintf('%1.0i',T)];
    savefile= '/home/veronika/synced/transfer_result/signal/sign/';
    save([savefile,savename],'xp_sign')
    clear all
    
end
        
