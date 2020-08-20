
clear all 
close all
%clc

place=1;

if place==1
    saveres=0;
    pltfig=1;
else
    saveres=1;
    pltfig=0;
end

nperm=1000;

tau=20;                                                                % choose between [10,20,50]
K=400;

%%
ba=1;
period=2;

namea={'V1','V4'};
namep={'tar','test'};
namebeh={'different','same'};
nameg={'SG','G','IG'};
ng=length(nameg);   % number of groups

task=['compute LDR remove info for layers ', namea{ba}, ' K = '  sprintf('%1.0i',K)];
display(task)

% exponential kernel for convolution
L=100;                                                                         % length of the kernel                                  
tau_vec=0:L;                                                                   % support
lambda=1/tau;                                                                  % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/transfer_mc/'); % for spikes                                                                       
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/'); % for weights
addpath('/home/veronika/synced/transfer_result/weight/tags/');       % for nb of neurons in the layer   

if place ==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results
      
loadname=['weight_bac_',namea{ba} ,'_',namep{period},'_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_',namep{period},'_', sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['ncell_lay_',namea{ba}];
load(loadname3)

%% reconstruct signal layers

nbses=length(w_all);
ncv=size(w_all{1},1);

x_all=zeros(nbses,ng,2,K);           % (sessions, groups, behavior,time)

tic
parfor sess = 1:nbses
    
    %%
    weights=w_all{sess};
    
    ncl=ncell_lay(sess,:);
    s=cumsum([0,ncl]);
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false);             % both conditions
    
    %% randomly permute the label "same" and "different" for spike trains of neurons of other sign
    
    x_partial=zeros(ng,nperm,2,K);
    for r = 1:ng
        tagged = s(r) + 1 : s(r+1);
        
        [x_cvpb] = reconstruct_partialb_fun(weights,spike_one,kernel,J,nperm,tagged);
        x_partial(r,:,:,:)=x_cvpb;                                                                  
        
    end
    
    x_all(sess,:,:,:)=squeeze(nanmean(x_partial,2)); % average across permutations
    
end
toc

%% plot

if pltfig==1
    
    blue=[0,0.48,0.74];
    col={blue,'g'};
    xvec=1:K;
    
    figure('units','centimeters','Position',[0,0,24,18])
    for r=1:ng   
        
        x_show=squeeze(x_all(:,r,:,:)); % group
        
        
        y1=squeeze(nanmean(x_show(:,1,:))-nanstd(x_show(:,1,:))./sqrt(nbses))';
        y2=squeeze(nanmean(x_show(:,1,:))+nanstd(x_show(:,1,:))./sqrt(nbses))';
        
        z1=squeeze(nanmean(x_show(:,2,:))-nanstd(x_show(:,2,:))./sqrt(nbses))';
        z2=squeeze(nanmean(x_show(:,2,:))+nanstd(x_show(:,2,:))./sqrt(nbses))';
        
        subplot(ng,1,r)
        hold on
        patch([xvec fliplr(xvec)], [y1 fliplr(y2)], col{1},'FaceAlpha',0.3,'EdgeColor',col{1})
        patch([xvec fliplr(xvec)], [z1 fliplr(z2)], col{2},'FaceAlpha',0.3,'EdgeColor',col{2})
        plot(zeros(K,1),'k')
        hold off
        
        box off
        axis([0,K,-0.004,0.004])
        text(1.02,0.5,nameg{r}, 'units','normalized')
        
        if r==1
            for ii=1:2
                text(0.05,0.8+(ii-1)*0.1,namebeh{ii},'units','normalized','color',col{ii})
            end
        end
        
        xlabel('time (ms)')
        ylabel('pop. signal')
    
    end
    
end

%% save

if saveres==1
    savename=['layer_remove_info_',sprintf('%1.0i',K)];
    savefile= '/home/veronika/synced/transfer_result/signal/layer/';
    save([savefile,savename],'x_all')
    
    clear variables
    
end
        
       
