% angle between population vectors

clear all 
close all
clc 

saveres=1;

Krange=[300,400,500];

task='compute angle between vectors w_S+C/w_C';
disp(task)

%% load results

addpath('/home/veronika/synced/transfer_result/weight/weight_comp2/')

loadname0='weight_comp_300';
load(loadname0)

nbses=size(w_comp,2);
ncv=size(w_comp{1},1);

%%

angle=zeros(3,nbses);

for t=1:3
    
    K=Krange(t);
    
    clear w_comp
    loadname=['weight_comp_',sprintf('%1.0i',K)];
    load(loadname)
    
    %% compute angle in cv
    
    for sess=1:nbses
        
        alpha=zeros(ncv,1);
        
        for cv=1:ncv
            
            x=w_comp{1,sess}(cv,:);
            y=w_comp{2,sess}(cv,:);
            alpha(cv)=acos(sum(x.*y));                                      % angle in radians
            
        end
        
        angle(t,sess)=mean(alpha);                                          % mean across cv
        
        
    end
    
end

alpham=mean(angle,2); % mean across sessions

%% angle permuted

nperm=1000;
anglep=zeros(3,nbses,nperm);

for t=1:3
    
    K=Krange(t);
    
    clear w_comp
    loadname=['weight_comp_',sprintf('%1.0i',K)];
    load(loadname)
    
    %% compute angle in cv
    
    for sess=1:nbses
        
        alphap=zeros(ncv,nperm);
        N=size(w_comp{1,1},2);
        for cv=1:ncv
            
            x=w_comp{1,sess}(cv,:);
            y=w_comp{2,sess}(cv,:);
            
            a=max(abs(x));
            b=max(abs(y));
            
            for perm=1:nperm
                
                xp =- a + (a + a).*rand(N,1);
                yp =- b + (b + b).*rand(N,1);
                
                alphap(cv,perm)=real(acos(sum(xp.*yp)));                          % angle in radians
            end
            
        end
        
        anglep(t,sess,:)=mean(alphap);                                      % mean across cv
        
        
    end
    
end

alpha_mp=squeeze(mean(anglep,2)); % mean across sessions

%% permutation test

pval=zeros(3,1);
for t=1:3
    
    x=alpham(t);
    x0=alpha_mp(t,:);
    pval(t)=sum(x>x0)/nperm;
    
end

%%
if saveres==1
    
    savename='w_angle';
    savefile='/home/veronika/synced/transfer_result/weight/w_angle/';
    save([savefile,savename],'alpham','alpha_mp','pval','Krange')
    %clear all
    
end




