function [wtilde] = weights_label_fun(s1,s2,nfold,Cvec,y_train)

format long
%%%%%%%%%%%%%%%%% computes linear svm on spike counts 
% details:
% binary classification of high-dimensional activity profiles of simultaneously active neurons; input dimensionality: N-by-ntrial
% all trials are used
% standardized the data
% estimates the best C-parameter with 10-fold cross-validation
% computes feature weights from support vectors

%% standardize

N=size(s1,2);

mat=cat(1,s1,s2);
stds=repmat(std(mat),size(mat,1),1);
means=repmat(mean(mat),size(mat,1),1); % regularize spike counts by subtracting the mean and dividing by std, for each neuron
mat_train=(mat-means)./stds;

%% select C-parameter

rp=randperm(size(mat_train,1));
matnp=mat_train(rp,:);
y_train=y_train(rp);

Ntr=floor(size(matnp,1)/nfold);         % nb trials in a split
bac_c=zeros(length(Cvec),nfold);

for j=1:length(Cvec)                    % cross-validation for splits into training and test data
    C=Cvec(j);
    
    for m=1:nfold                       % n-fold cross-validation for splitting again the training data into training and validation set
        
        xtest=matnp(1+(m-1)*Ntr:m*Ntr,:);                                       % data for testing
        xtrain=[matnp(1:(m-1)*Ntr,:);matnp(m*Ntr+1:end,:)];                     % data for training
        labc_train=[y_train(1:(m-1)*Ntr);y_train(m*Ntr+1:end)];                 % label training
        labc_test=y_train(1+(m-1)*Ntr:m*Ntr);                                   % label testing
        
        try
            
            svmstruct=svmtrain(xtrain,labc_train','kernel_function','linear','boxconstraint',C,'autoscale',false); % train linear svm
            class=svmclassify(svmstruct,xtest);                                 % classify the validation data
            
            tp =length(find(labc_test==1 & class==1));                          % TruePos
            tn =length(find(labc_test==-1 & class==-1));                        % TrueNeg
            fp =length(find(labc_test==-1 & class==1));                         % FalsePos
            fn =length(find(labc_test==1 & class==-1));                         % FalseNeg
            
            % balanced accuracy
            if tp==0
                bac_c(j,m)=tn./(tn+fp);
            elseif tn==0
                bac_c(j,m)=tp./(tp+fn);
            else
                bac_c(j,m) =((tp./(tp+fn))+(tn./(tn+fp)))./2;
            end
        catch
            bac_c(j,m)=0;
            
        end
    end
    
end

[~,idx]=max(mean(bac_c,2));
C=Cvec(idx);


%% Extract feature weights

try  
    %% train the SVM classifier with selected C-parameter  
    
    svmstruct=svmtrain(matnp,y_train,'kernel_function','linear','boxconstraint',C, 'autoscale',false);   % 
                                                                 
    %% compute weights
    
    svt=svmstruct.SupportVectorIndices;                                                   % indices of used trials
    y=y_train(svt);                                                                       % label of used trials
    
    wi=zeros(length(svmstruct.Alpha),N);                                                  % svm.SupportVectors~(n_used_trials,n_neurons)
    for k=1:length(svmstruct.Alpha)                                                       % svm.Alpha~ Lagrange multiplier
        wi(k,:)=svmstruct.Alpha(k,1)*y(k,1)*svmstruct.SupportVectors(k,:);                % length(svm.Alpha) ~ number of support vectors
    end
    
    w=sum(wi,1)';  % feature weights of the SVM                  
    wtilde=w./sqrt(sum(w.^2));
                                                        % svm.Bias~Intercept of the hyperplane that separates the two groups in the normalized data space
catch
    
    wtilde=zeros(1,N);
    disp('No convergence!')
    
end


