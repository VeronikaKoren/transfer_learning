function [x_cvpb] = reconstruct_partialb_fun(weights,spike_one,kernel,J,nperm,tagged)

%% Low-dimensional reconstruction of spike trains

% computes low-dimensional reconstruction of spike trains as the weighted
% sum of spikes, convolved with an exponential kernel

% computes LDR in trials and in cross-validations and averages across
% trials and across cvs

%%

format short
ncv=size(weights,1);                                                        % number of neurons
K=size(spike_one{1},3);                                                     % number timesteps

%% compute LDR

xcv=zeros(ncv,nperm,2,K);
idx_keep=tagged;

for cv=1:ncv
    
    spikone=spike_one{cv};
    w=weights(cv,:);
    xsp=zeros(nperm,2,K);
    for perm=1:nperm
        
        %%%%%%%%%%%%% permute trials
        
        spikes_perm=spikone(randperm(sum(J)),:,:);
        spikes_perm(:,idx_keep,:)=spikone(:,idx_keep,:);                                                      % neurons with tag have regular spike counts
        spikes=cell(2,1);
        spikes{1}=spikes_perm(1:J(1),:,:);
        spikes{2}=spikes_perm(J(1)+1:end,:,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x_full=cell(2,1);
        for c=1:2
            
            o=single(spikes{c});
            xsig=zeros(J(c),K);
            for trial=1:J(c)
                
                oj=squeeze(o(trial,:,:));                                                       % reconstruct
                xsig(trial,:)=conv(w*oj,kernel./sum(kernel),'same');
                
            end
            x_full{c}=xsig;
        end
        
        % subtract the running mean
        z=mean(cat(1,x_full{1},x_full{2}));
        repz=cellfun(@(x,y) repmat(x,y,1),[{z};{z}],[{J(1)};{J(2)}],'UniformOutput',false);
        nsig=cellfun(@(x,z) x-z,x_full,repz,'UniformOutput',false);                             % deviation from the mean
        msig=cellfun(@mean, nsig,'UniformOutput',false);                                        % average across trials
        
        %%%%%%%%%%%%%%%%%%%
        xsp(perm,1,:)=msig{1};
        xsp(perm,2,:)=msig{2};
    end
    xcv(cv,:,:,:)=xsp;
    
    
end
%% average across cv

x_cvpb=squeeze(nanmean(xcv));

end

