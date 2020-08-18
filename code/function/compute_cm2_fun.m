function [r_cm] = compute_cm2_fun(spike_train,D,u)
%compute cross-mapping for pairs of neurons in trials

                                                                   
J=size(spike_train,1);                                                      % number of trials    
N=size(spike_train,2);                                                      % number of neurons
K=size(spike_train,3);                                                      % number of time steps (length of the signal)
Np=(N^2-N)/2;                                                               % number of pairs

E=D+1;                                                                      % number of nearest neighbors
L=K-D+1;                                                                    % clipped length of the signal (length - delay)                                              

%% pre-compute pairs of spike trains

sts=cell(Np,2);

counter=0;
for n=1:N-1
    
    st1=squeeze(spike_train(:,n,:));
    
    for m=n+1 :N
        
        counter=counter+1;
        
        st2=squeeze(spike_train(:,m,:));
        sts{counter,1}=st1;
        sts{counter,2}=st2;
        
    end
end

%% compute cross-mapping

r_cm=zeros(2,Np);

for pid=1:Np
    
    st1=sts{pid,1};
    st2=sts{pid,2};
    
    ryx=zeros(J,2);
    
    for j = 1:J
        
        x0=conv(st1(j,:),u,'same');
        y0=conv(st2(j,:),u,'same');
        
        %%
        
        X=zeros(D,L);
        Y=zeros(D,L);
        
        X(1,:)=x0(D:end);
        Y(1,:)=y0(D:end);
        
        for d=2:D
            X(d,:)=x0((D:end)-(d-1));
            Y(d,:)=y0((D:end)-(d-1));
        end
        
        
        %% random weighting of coordinates
        rw=randn(D);                                                                % random weights (normal distribution?)
        
        Xrw=zeros(D,L);
        Yrw=zeros(D,L);
        for t=1:L
            Xrw(:,t)=rw'*X(:,t);
            Yrw(:,t)=rw'*Y(:,t);
        end
        
        %% matrix of similarity
        
        sim_x=zeros(L,L);
        sim_y=zeros(L,L);
        
        for t1=1:L
            for t2=1:L
                
                % neuron 1
                distx=sqrt(sum((Xrw(:,t2)-Xrw(:,t1)).^2));                              % distance of points in random delay coordinates
                nrm=sqrt(sum((Xrw(:,1)-Xrw(:,t1)).^2));                                 % normalization
                sim_x(t1,t2)=exp(-distx/nrm);                                           % similarity of points in random delay coordinates
                
                % neuron 2
                disty=sqrt(sum((Yrw(:,t2)-Yrw(:,t1)).^2));
                nrmy=sqrt(sum((Yrw(:,1)-Yrw(:,t1)).^2));
                sim_y(t1,t2)=exp(-disty/nrmy);
            end
        end
        
        %sim_x=sim_x-eye(size(sim_x));                                                  % remove the diagonal
        %sim_y=sim_y-eye(size(sim_y));
        
        %% compute nearest neighbors
        
        nnx=zeros(E,L);
        nny=zeros(E,L);
        
        for t=1:L
         
            [~,idx]=sort(sim_x(t,:));                                             % indexes of E nearest neighbors from the similarity matrix; sorts indexes from smallest to largest
            %nnx(t,:)=idx(end-(E-1):end)                                          % from 1st to E NN [training]
            nnx(:,t)=idx(end-E:(end-1));                                         % from 2nd to E+1 NN [testing]
            
            [~,idy]=sort(sim_y(t,:));
            %nny(t,:)=idy(end-(E-1):end);
            nny(:,t)=idy(end-E:(end-1));
            
        end
        
        %% reconstruction using nearest neighbors of the other signal
        
        crossy=zeros(E,D,L); % (neighbors, dimension, timesteps)
        crossx=zeros(E,D,L);
        for e=1:E
            
            crossy(e,:,:)=Yrw(:,nnx(e,:));                                             % cross-mapping of Y using neighbors from X
            crossx(e,:,:)=Xrw(:,nny(e,:));                                             % cross-mapping of X using neighbors from Y
            
        end
        
        Yhat=squeeze(mean(crossy,1));                                                   %mean accross E nearest neighbors
        Xhat=squeeze(mean(crossx,1));
        
        %% linear correlation between the original projection and the cross-projection
        
        ryx(j,1)=corr(Yrw(1,:)',Yhat(1,:)');                                            % collect across trials
        ryx(j,2)=corr(Xrw(1,:)',Xhat(1,:)');
        
    end
    
    r_cm(:,pid)=nanmean(ryx);                                                        % average across trials and collect across pairs
    
end


end

