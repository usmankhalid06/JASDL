function [H,Z,D,X]= JASDL_sim(Y,Hq,Zq,Dq,Xq,nIter,K,W1,W2,zeta,lambda,alpha,fac,sub)
    A = eye(size(Dq,2),K); 
    S = eye(K,size(Hq,1)); 
    B = zeros(K,size(Xq,1));
    T = zeros(size(Zq,2),K);
    D = Dq*A;
    
%     fprintf('Iteration:     ');    
    for iter=1:nIter
%         fprintf('\b\b\b\b\b%5i',iter);   
        D = Dq*A;
        X = B*Xq;
        H = S*Hq;
        Z = Zq*T;   

        if sub<=3
            for jjj=1:fac(1)
                [~,ind4(jjj)] = max(abs(corr(Dq(:,jjj),H')));
            end
        else
            for jjj=[1:fac(1) fac(2)]
                [~,ind4(jjj)] = max(abs(corr(Dq(:,jjj),H')));
            end            
        end
        
        
        for j=1:size(D,2)
            X(j,:) = 0; A(:,j) = 0; B(j,:) = 0; S(j,:) = 0; T(:,j) = 0; 
            E = Y-D*X;
            if any(j==ind4)
                zeta2 = W1(1);  zeta3 = W2(1);
            else
                zeta2 = W1(2);  zeta3 = W2(2);
            end
            

            xk = alpha*D(:,j)'*E+(1-alpha)*H(j,:)*E; 
            thr = lambda./abs(xk);
            xkk = sign(xk).*max(0, bsxfun(@minus,abs(xk),thr/2));   
            [~,bb1]= sort(abs(Xq*xkk'),'descend');
            ind1 = bb1(1:zeta3);
            B(j,ind1)= xkk/Xq(ind1,:); %xkk*Xq(ind1,:)'/(Xq(ind1,:)*Xq(ind1,:)');
            X(j,:) = B(j,:)*Xq;
            
            xk2 = E'*D(:,j);
            thr2 = lambda./abs(xk2);
            xkk2 = sign(xk2).*max(0, bsxfun(@minus,abs(xk2),thr2/(2*(1-alpha))));
            [~,bb1]= sort(abs(Zq'*xkk2),'descend');
            ind1 = bb1(1:zeta);
            T(ind1,j)= Zq(:,ind1)\xkk2;   %xkk*Zq(ind1,:)'/(Zq(ind1,:)*Zq(ind1,:)');
            Z(:,j) = Zq*T(:,j);
            
            rInd1 = find(X(j,:)); 
            if (length(rInd1)<1)
                [~,indd]= max(sum(Y-Dq*A*X.^2));
                D(:,j)= Y(:,indd)/norm(Y(:,indd));
            else
                tmp1 = X(j,rInd1)*E(:,rInd1)';
                [~,bb3]= sort(abs(Hq*tmp1'),'descend');
                ind3 = bb3(1:zeta);
                S(j,ind3)= tmp1/Hq(ind3,:);      %((Hq(ind3,:)'*Hq(ind3,:))\Hq(ind3,:)')'*tmp4';   
                S(j,:) = S(j,:)./norm(S(j,:)*Hq);   
                H(j,:) = S(j,:)*Hq;
                
                
                mag = X(j,rInd1)*X(j,rInd1)';
                tmp2 = (alpha*E(:,rInd1)*X(j,rInd1)'+(1-alpha)*E(:,rInd1)*Z(rInd1,j))/(1-alpha+(alpha*mag));  
                [~,bb2]= sort(abs(Dq'*tmp2),'descend');
                ind2 = bb2(1:zeta2);
                A(ind2,j)= Dq(:,ind2)\tmp2;     %(Dq(:,ind2)'*Dq(:,ind2))\Dq(:,ind2)'*tmp3;      
                A(:,j) = A(:,j)./norm(Dq*A(:,j));    
                D(:,j) = Dq*A(:,j);
            end
        end
        
    end
    H =H'; Z=Z';
end



