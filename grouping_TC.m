function [tmpDc,tmpXc] = grouping_TC(SM,T,S,nS)

for k =1:size(SM,1)
    tmpD=[];
    tmpX=[];
    for j=1:nS
        [~,icc] = max(abs(corr(abs(SM(k,:)'),abs(S(:,:,j)')))); %max(abs(corr(abs(SM(j,:)'),abs(Zs'))));
        tmpD = [tmpD T(:,icc,j)];
        tmpX = [tmpX;S(icc,:,j)];
    end
    [tmpDc(:,k),tt1,tt2] = svds((tmpD*tmpX)/nS,1);
    tmpDc(:,k) = tmpDc(:,k)/norm(tmpDc(:,k));
    tmptmpXc = tt1*tt2';
    tmpXc(k,:) = sign(tmptmpXc).*max(0, bsxfun(@minus,abs(tmptmpXc),0.01/2));
end