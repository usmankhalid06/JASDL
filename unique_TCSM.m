function [tmpD,tmpX] = unique_TCSM(tmpTC,Zt,Zs,M)

for j=1:M
    [~,icc] = max(abs(corr(tmpTC(:,j),Zt(:,:,j))));
    tmpD(:,j) = Zt(:,icc,j);
    tmpX(j,:) = Zs(icc,:,j);
end
