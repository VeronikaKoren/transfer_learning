% reconstruct signal with bursty neurons

clear all 
close all
clc

place=1;

if place==1
    saveres=0;
    pltfig=1;
    
else
    pltfig=0;
    saveres=1;
end

nperm=1000;

%%
ba=1;
period=2;

tau=20;                                                                % choose between [10,20,50]
K=400;

%%

namebeh={'different','same'};
namea={'V1','V4'};
nameb={'bursty','nonbursty'};
namep={'tar','test'};

task=['compute LDR for bursty and nonbursty neurons ', namea{ba}, ' K = '  sprintf('%1.0i',K)];
display(task)

% exponential kernel for convolution
L=100;                                                                                % length of the kernel                                  
tau_vec=0:L;                                                                              % support
lambda=1/tau;                                                                   % time constant
kernel=exp(-lambda.*tau_vec); 

%% add path

addpath('/home/veronika/synced/transfer_result/input/transfer_mc/') 
addpath('/home/veronika/synced/transfer_result/weight/tags/')
addpath('/home/veronika/synced/transfer_result/weight/weight_bac/')

if place ==1
    addpath('/home/veronika/Dropbox/transfer/code/function/')
else
    addpath('/home/veronika/transfer/code/function/')
end

%% load results

loadname=['weight_bac_',namea{ba},'_' ,namep{period},'_', sprintf('%1.0i',K)];
load(loadname,'w_all')

loadname2=['input_',namea{ba},'_', namep{period},'_',sprintf('%1.0i',K)];
load(loadname2,'spikes_inm','spikes_cnm')

loadname3=['tag_bursty_', sprintf('%1.0i',K)];
load(loadname3,'bursty')

%% reconstruct signal for bursty and non-bursty neurons

nbses=size(spikes_cnm,1);
ncv=size(spikes_cnm,2);

x1=zeros(nbses,nperm,2,K);
x2=zeros(nbses,nperm,2,K);

tag=[1,0];   % bursty, non-bursty

tic
parfor sess = 1:nbses
    
    %%
    weights=w_all{sess};
    
    spike_diff=(spikes_cnm(sess,:)');
    spike_same= repmat(spikes_inm(sess),ncv,1);
    J=[size(spike_diff{1,1},1), size(spike_same{1,1},1)];
    
    spike_one=cellfun(@(x,y) cat(1,x,y), spike_diff,spike_same,'UniformOutput', false); % both conditions
    
    %% randomly permute the label "same" and "different" for spike trains of neurons of the other group
    
    for r = 1:2
        
        tagged=find(bursty{sess}==tag(r));                                             % idx to keep
        
        x_perm=zeros(nperm,2,K);
        for perm=1:nperm
            
            rp=randperm(sum(J));                                                        % random order of trials
            spikes_perm=cellfun(@(x) x(rp,:,:), spike_one,'UniformOutput', false);      % random order to all spike trains    
           
            for cv=1:ncv                                                                % give correct order of trials to the group
                spikes_perm{cv}(:,tagged,:)=spike_one{cv}(:,tagged,:);
            end
             
            spikes=cell(ncv,2);
            spikes(:,1)=cellfun(@(x) x(1:J(1),:,:), spikes_perm, 'UniformOutput', false);               % put in form cell{ncv,2}
            spikes(:,2)=cellfun(@(x) x(J(1)+1:end,:,:), spikes_perm, 'UniformOutput', false);
            
            %%
            
            [x_cv] = reconstruct_fun(weights,spikes,kernel);
            x_perm(perm,:,:)=x_cv;                                      % collect across permutations
            
        end
        
        if r==1
            x1(sess,:,:,:)=x_perm;                                      % collect results across sessions
        else
            x2(sess,:,:,:)=x_perm;
        end
         
    end
end
toc

%%

x_bursty=cell(2,1);   % put in a cell
x_bursty{1}=x1; 
x_bursty{2}=x2;

%% plot

if pltfig==1
    
    x=1:K;
    blue=[0,0.48,0.74];
    col={blue,'g'};
    
    figure('units','centimeters','Position',[0,0,24,18])
    for r=1:2
        
        x_show=squeeze(nanmean(x_bursty{r},2)); % average across permutations
        
        subplot(2,1,r)
        hold on
        for ii=1:2
           y1=squeeze(nanmean(x_show(:,ii,:)) - nanstd(x_show(:,ii,:))./sqrt(nbses))';
           y2=squeeze(nanmean(x_show(:,ii,:)) + nanstd(x_show(:,ii,:))./sqrt(nbses))';  
           patch([x fliplr(x)], [y1 fliplr(y2)], col{ii},'FaceAlpha',0.3,'EdgeColor',col{ii})
        end
        
        plot(zeros(K,1),'k')
        hold off
        box off
        axis([0,K,-0.004,0.004])
        text(1.05,0.5,nameb{r}, 'units','normalized')
        
        if r==1
            text(0.1,0.95,namebeh{1},'units','normalized','color',col{1})
            text(0.1,0.87,namebeh{2},'units','normalized','color',col{2})
        else
            xlabel('time (ms)')
        end
        ylabel('pop. signal')
    
    end
    
end

%% save

if saveres==1
   
    savename=['bursty_ri_',sprintf('%1.0i',K)];
    savefile= '/home/veronika/synced/transfer_result/signal/bursty/';
    save([savefile,savename],'x_bursty')
    
    %clear all
    
end
        
       
