%% power spectrum of spike trains
% with permuted spike timing to get the null model

close all
clear all
clc

saveres=0;

nperm=1000;
K=400;  % length of the time window

%%
ba=1;
period=2;

namea={'V1','V4'};
namep={'target', 'test'};
namei={'s+c','c'};

start_vec=[200,500];
start=start_vec(period);                                        % start time

%% PS parameters

dt = 0.001;     % in seconds
fs = 1/dt;      % sampling frequency in Hz
Nfreq=fs/2;     % Nyquist frequency=half of the sampling rate at which the signal is sampled.
NW=3;           % time-bandwidth product (uses 2*NW-1 tapers)

n=0:K-1;
w=0.5*(1-cos((2*pi).*n./K));    % Hann window


nfft= floor(K/2)+1;             % number of fft
delta_f= Nfreq./(nfft-1);       % resolution
display(delta_f,'resolution in freq')

%% get spike trains

addpath '/home/veronika/synced/transfer_result/input/spike_train/';

loadname=['spike_train_s+c_',namea{ba},'_',namep{period}];
load(loadname);

stime=cellfun(@(x) x(:,:,start:start+K-1),spiketrain,'UniformOutput', false);
strain=cellfun(@(x,y) cat(1,x,y),stime(:,1),stime(:,2),'UniformOutput',false); % concatenate conditions

nbses=size(spiketrain,1);

%% compute PS

task=['compute power spectrum with permuted timing in ' namea{ba} namep{period}];
disp(task)

tic
psp=cell(nbses,1);

for sess=1:nbses
    
    spike_sess=single(strain{sess});
    
    N=size(spike_sess,2);
    J=size(spike_sess,1);
    frate=mean(mean(spike_sess,1),3).*1000;             % mean firing rate
    
    %% compute PS on permuted spike timing (PS method with Slepian tapers)
    
    ps_perm=zeros(nperm,N,floor(nfft/2)+1);
    for perm=1:nperm
        
        x_perm=spike_sess(:,:,randperm(K));                % permute spike timing
        
        ps_all=zeros(N,floor(nfft/2)+1);
        for ii=1:N
            
            spikes=squeeze(x_perm(:,ii,:));
            
            pxx=pmtm(spikes',NW,nfft,fs);                  % compute PS in every trial
            ps=(pxx./(dt^2*K/nfft))';                      % normalization
            ps_all(ii,:)=mean(ps)./frate(ii);              % average across trials and normalize with the firing rate
            
        end
        
        ps_perm(perm,:,:)=ps_all;                          % collect across permutations
        
    end
    %%
    
    psp{sess}=ps_perm;                          % collect across sessions
    
end
toc

f=0:delta_f*2:fs/2;

%%

if saveres==1
    savename=['ps_poisson_',namea{ba},namep{period},'_',sprintf('%1.0i',K) ];
    savefile='/home/veronika/synced/struct_result/basic_stat/ps/';
    save([savefile,savename],'psp','T','f')
end

