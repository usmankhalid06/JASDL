clear all
close all
clc




%% loading SimTB sources and defining other parameters. To generate a different dataset you
%% can download the SimTB on your computer

load sources_SimTB
% TC_gw is the ground truth temporal sources defined by the first subject and common sources 
% SM_gw in the ground truth spatial maps defined by the first subject and common sources
% iTC_sw contains TC_sw in 3D arrays
% iSM_sw contains SM_sw in 3D arrays

N = size(TC_gw,1);
V = size(SM_gw,2);
nCom = 6; % common components 
K = 12; % total number of components to estimate
srcs = 7; % srcs per subject 6 common +1 unique 
nIter = 30; % number of algorithm iterations
Dp = odctdict(N,N); % generating DCT bases
Dp = Dp(:,2:end); % removing the constant
tstd  = sqrt(1.2); %% noise vector generation for time sources
sstd  = sqrt(0.02); %0.0025
M = 6; %6 subjects
nA = 4; %number of participating algorithms
 
%% Generating the data using spatial and temporal sources
vec = 6; 
rng('default'); 
% rng('shuffle') 
rng(50,'philox') 
for sub=1:M
    TC_sw{sub} = reshape(iTC_sw(sub,:,:),N,srcs); %TC_sw are the subject-wise temporal sources for all 6 subjects 
    SM_sw{sub} = reshape(iSM_sw(sub,:,:),srcs,V); % SM_sw are the subject-wise spatial sources for all 6 subjects
    Y{sub} = (TC_sw{sub}+tstd*randn(N,nCom+1))*(SM_sw{sub}+sstd*randn(nCom+1,V));
    Y{sub} = Y{sub}-repmat(mean(Y{sub}),size(Y{sub},1),1);
end


%% Applying algorithms on the generated data
%% ICA algorithm
[Zt{1},Zs{1},D{1},X{1}] = cgICA_sim(Y,M,K,K);

%% ssBSS_pre method
params1.K = 14; %14
params1.P = 7;
params1.lam1 = 3;
params1.zeta1 = 60;
params1.Kp = 150;
params1.nIter = nIter;
params1.alpha = 10^-8;
for i=1:M
    [Zt{2}(:,:,i),Zs{2}(:,:,i),~,~]=ssBSS_pre(Y{i},Dp,params1,TC_sw{i},SM_sw{i}); %_mod
end

%% JASDL 
% concatenating ssBSS spatial and temporal components
TT = [Zt{2}(:,:,1) Zt{2}(:,:,2) Zt{2}(:,:,3) Zt{2}(:,:,4) Zt{2}(:,:,5) Zt{2}(:,:,6)];
SS = [Zs{2}(:,:,1);Zs{2}(:,:,2);Zs{2}(:,:,3);Zs{2}(:,:,4);Zs{2}(:,:,5);Zs{2}(:,:,6)];
SS = abs(SS);
for jjj=1:size(SS,1)
    SS(jjj,:) =(SS(jjj,:) - min(SS(jjj,:))) / (max(SS(jjj,:)) - min(SS(jjj,:)));
end      

% using biological priors in case of real fMRI data, but in synthetic we using source TCs and SMs
TC_gw2=TC_gw*diag(1./sqrt(sum(TC_gw.*TC_gw))); 
SM_gw2 = abs(SM_gw);
for jjj=1:size(SM_gw,1)
    SM_gw2(jjj,:) =(SM_gw2(jjj,:) - min(SM_gw2(jjj,:))) / (max(SM_gw2(jjj,:)) - min(SM_gw2(jjj,:)));
end

% preparing spline bases
PSF_sim3 = spcol([0, 0, 0, linspace(0,1,34), 1, 1, 1], 4, linspace(0,1,50))';  
PSF_sim4 = kron(PSF_sim3, PSF_sim3); 
for jjj=1:size(PSF_sim4,1)
    PSF_sim4(jjj,:) =(PSF_sim4(jjj,:) - min(PSF_sim4(jjj,:))) / (max(PSF_sim4(jjj,:)) - min(PSF_sim4(jjj,:)));
end

alpha = 0.45; %0.3
W1 = [6,30]; %choosing different sparsity for known temporal dynamics and different for unknown (paper's equation 5 Omega1)
W2 = [6,120]; %choosing different sparsity for known spatial dynamics and different for unknown (paper's equation 5 Omega2)
zeta = 30; % Sparsity parameter for zeta in paper's equation 5
lambda = 6; % Saprsity parameter for soft thresholding\
% parameter F is different here than real data because we are assuming that
% all sources are known for last 3 subjects where as the unique source is
% unknown for the first 3 subjects
for i =1:M
    Hq = TT';  Zq = SS'; %prepared using ssBSS method  
    Dq2 = [TC_gw2(:,[1:6 10:12]) Dp(:,1:150)];  %concatenating MHRs with DCT bases
    Xq2 = [SM_gw2([1:6 10:12],:); PSF_sim4];  %concatenating FBNs/IBNs with Splines
    if i<=3
        F = [6 6];
        [Zt{3}(:,:,i),Zs{3}(:,:,i),Zt{4}(:,:,i),Zs{4}(:,:,i)...
            ]= JASDL_sim(Y{i},Hq,Zq,Dq2,Xq2,nIter,K,W1,W2,zeta,lambda,alpha,F,i);   
    else
        F = [7 7];
        [Zt{3}(:,:,i),Zs{3}(:,:,i),Zt{4}(:,:,i),Zs{4}(:,:,i)...
            ]= JASDL_sim(Y{i},Hq,Zq,Dq2,Xq2,nIter,K,W1,W2,zeta,lambda,alpha,F,i); 
    end
end

         

%% Analysis
for j =2:nA
    [D{j}(:,1:6),X{j}(1:6,:)] = grouping_TC(SM_gw(1:6,:),Zt{j},Zs{j},M);
    [D{j}(:,7:12),X{j}(7:12,:)] = unique_TCSM(TC_gw(:,7:12),Zt{j},Zs{j},M);
end

for a =1:nA
    [sD{a+1},sX{a+1},ind]=sort_TSandSM_spatial(TC_gw,SM_gw,D{a},X{a},K);
    for i =1:nCom+6
        TCcorr_gw(i,a+1) =abs(corr(TC_gw(:,i),D{a}(:,ind(i))));
        SMcorr_gw(i,a+1) =abs(corr(SM_gw(i,:)',X{a}(ind(i),:)'));
    end
end
            

sD{1} = TC_gw; 
sX{1} = SM_gw;
tt1 = TCcorr_gw;
tt2 = SMcorr_gw;
figure; my_subplots_github(nA+1,K,sqrt(V),sqrt(V),tt1',tt2',sD,sX);
 

TP = zeros(nA,12); FP = zeros(nA,12); FN = zeros(nA,12); F_score = 0; thr = 0.02;
for i =1:nA
    for jjj=1:12
        SM_gw4(jjj,:) =SM_gw(jjj,:)/norm(SM_gw(jjj,:)); 
        [~, indd(jjj)]  = max(abs(corr(abs(SM_gw4(jjj,:)'),abs(X{i}'))));
        XX{i}(indd(jjj),:) = X{i}(indd(jjj),:)/norm(X{i}(indd(jjj),:));
        TP(i,jjj) = TP(i,jjj) +sum(abs(SM_gw4(jjj,:))>=thr & abs(XX{i}(indd(jjj),:))>=thr);
        FP(i,jjj) = FP(i,jjj) +sum(abs(SM_gw4(jjj,:))<=thr & abs(XX{i}(indd(jjj),:))>=thr);
        FN(i,jjj) = FN(i,jjj) +sum(abs(SM_gw4(jjj,:))>=thr & abs(XX{i}(indd(jjj),:))<=thr);   
    end
    F_score(i) = (2*sum(TP(i,:)))/(2*sum(TP(i,:))+sum(FP(i,:))+sum(FN(i,:)));
end


TC_corr = mean(tt1);
SM_corr = mean(tt2);
F_score =F_score;    



%% Printing results
% Print the header
fprintf('Following are the results for ICA, ssBSS, Proposed (H,Z), Proposed (D,X):\n');
fprintf('--------------------------------------------------------------------------------------------------\n');

% Table header
fprintf('%-60s', 'Metric');
fprintf('%8s   %8s   %8s %8s\n', 'ICA', 'ssBSS', 'Proposed(H,Z)', 'Proposed(D,X)');

% Correlation between source TCs and recovered TCs
fprintf('%-60s', 'Correlation between source TCs and recovered TCs:');
fprintf('  %8.4f  %8.4f   %8.4f      %8.4f\n', TC_corr(2:end));

% Correlation between source SMs and recovered SMs
fprintf('%-60s', 'Correlation between source SMs and recovered SMs:');
fprintf('  %8.4f  %8.4f   %8.4f      %8.4f\n', SM_corr(2:end));

% Fscores for recovered SMs
fprintf('%-60s', 'Fscores for recovered SMs:');
fprintf('  %8.4f  %8.4f   %8.4f      %8.4f\n', F_score(1:end));

fprintf('--------------------------------------------------------------------------------------------------\n');