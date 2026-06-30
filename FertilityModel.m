function [x,beta,betaA] = FertilityCoef_P

load(char(pATh + "F.mat"),'x','F');
Z               = 0; %This is to incorporate additional covariates, but not for now. 
y               = log(F./sum(F,1)); %To have regular age patterns of fertility and estimate a model in logs. The log transformation applies to mortality data either using nMx or q(x).
[~,s,~]         = svd(y*y');
R               = sum(cumsum(diag(s))/sum(diag(s)) < (1 - 10^-5)); %To reduce the data, becuase rates are correlated (see the note at the end of the word file). 

yA              = y(:,2:2:end); %To represent the secondary dataset, i.e., DHS inputs in your case.
[u,s,v]         = svds(yA,R);
UA              = u;
YA              = (v*s)';             
XA              = ones(1,size(YA,2));
YYA             = YA*YA';
XXA             = XA*XA';
XYA             = XA*YA';
thetaA          = XXA\XYA;
EEA             = YYA - XYA'*thetaA - thetaA'*XYA + thetaA'*XXA*thetaA;
[u,~,~]         = svds(UA*EEA*UA',R);
UolsA           = u;

y               = y(:,1:2:end);  %To represent be the prior dataset, i.e., HMD inputs in your case.
[u,s,v]         = svds(y,R);
U               = u;
Y               = (v*s)';
X               = ones(1,size(Y,2));
YY              = Y*Y';
XX              = X*X';
XY              = X*Y';
theta           = XX\XY;
EE              = YY - XY'*theta - theta'*XY + theta'*XX*theta;
[u,~,~]         = svds(U*EE*U',R);
Uols            = u;

[ages,tables]   = size(Y);
[~,tablesA]     = size(YA);
I               = sparse(eye(ages));
priorB          = {zeros(1 + Z,ages),ones(1 + Z,ages)*10^-4};
priorPsi        = {ones(ages,1)*10^-2,ones(ages,1)*10^-2};
priorT          = {triu(NaN(ages)) + tril(zeros(ages),-1),triu(NaN(ages)) + tril(ones(ages),-1)*10^-2};
G               = sum(~isnan(priorT{1}));
G               = G(G > 0);
lisT            = find(~isnan(priorT{2}));

C               = sparse(zeros(ages^2,numel(lisT)));
for i = 1:numel(lisT)
    C(lisT(i),i) = 1;
end
B               = theta;
Psi             = sparse(diag(gamrnd(priorPsi{1},priorPsi{2})));
T               = eye(ages);
BA              = thetaA;
PsiA            = sparse(diag(gamrnd(priorPsi{1},priorPsi{2})));
TA              = eye(ages);

IterAT          = 5000;
BurnIn          = 250;
for i = 1:BurnIn + IterAT
    clc;
    i/(BurnIn + IterAT)
    EE           = YY - XY'*B - B'*XY + B'*XX*B;
    EEA          = YYA - XYA'*BA - BA'*XYA + BA'*XXA*BA;

    for j = 1:ages
        Psi(j,j)   = gamrnd(priorPsi{1}(j) + (tables)/2,1/(priorPsi{2}(j) + 1/2*T(:,j)'*EE*T(:,j)));
        PsiA(j,j)  = gamrnd(priorPsi{1}(j) + (tables)/2 + (tablesA)/2,1/(priorPsi{2}(j) + 1/2*T(:,j)'*EE*T(:,j) + 1/2*TA(:,j)'*EEA*TA(:,j)));
    end        
    
    L            = mat2cell(C'*kron(Psi,EE)*C + sparse(diag(priorT{2}(lisT))),G,G);
    V            = mat2cell(zeros(size(C,2)),G,G);
    for j = 1:numel(G)
        V{j,j} = sparse(pinv(full(L{j,j})));
        L{j,j} = chol(L{j,j})';        
    end
    L            = cell2mat(L);
    V            = cell2mat(V);
    Mu           = -C'*kron(Psi,EE)*reshape(I,[],1) + sparse(diag(priorT{2}(lisT))*priorT{1}(lisT));
    T(lisT)      = V*Mu + L'\normrnd(0,1,numel(lisT),1);
    Phi          = T*Psi*T';
    clear L V Mu

    L            = mat2cell(C'*(kron(Psi,EE) + kron(PsiA,EEA))*C + sparse(diag(priorT{2}(lisT))),G,G);
    V            = mat2cell(zeros(size(C,2)),G,G);
    for j = 1:numel(G)
        V{j,j} = sparse(pinv(full(L{j,j})));
        L{j,j} = chol(L{j,j})';        
    end

    L            = cell2mat(L);
    V            = cell2mat(V);
    Mu           = -C'*(kron(Psi,EE) + kron(PsiA,EEA))*reshape(I,[],1) + sparse(diag(priorT{2}(lisT))*priorT{1}(lisT));
    T(lisT)      = V*Mu + L'\normrnd(0,1,numel(lisT),1);
    PhiA         = TA*PsiA*TA';
    clear L V Mu

    V            = kron(Phi,XX) + diag(reshape(priorB{2},[],1));
    L            = chol(sparse(V))';
    V            = pinv(V);
    Mu           = kron(Phi,XX)*reshape(theta,[],1) + diag(reshape(priorB{2},[],1))*reshape(priorB{1},[],1);
    B            = V*Mu + L'\normrnd(0,1,numel(B),1);
    B            = reshape(B',Z + 1,ages);
    clear L V Mu

    V            = kron(Phi,XX) + kron(PhiA,XXA) + diag(reshape(priorB{2},[],1));
    L            = chol(sparse(V))';
    V            = pinv(V);
    Mu           = kron(Phi,XX)*reshape(theta,[],1) + kron(PhiA,XXA)*reshape(thetaA,[],1) + diag(reshape(priorB{2},[],1))*reshape(priorB{1},[],1);
    BA           = V*Mu + L'\normrnd(0,1,numel(BA),1);
    BA           = reshape(BA',Z + 1,ages);
    clear L V Mu

    [u,~,~]      = svds(U*pinv(full(Phi))*U',R);
    u            = u.*sign(sum(Uols.*u));
    u            = [U*B' u];
    k            = (y - u(:,1))'*u(:,2:end);

    [uA,~,~]     = svds(UA*pinv(full(PhiA))*UA',R);
    uA           = uA.*sign(sum(UolsA.*uA));
    uA           = [UA*BA' uA];
    kA           = (yA - uA(:,1))'*uA(:,2:end);

    if i > BurnIn
        for j = 1:size(u,2)
            beta{j}(:,i - BurnIn)  = u(:,j);
            betaA{j}(:,i - BurnIn) = uA(:,j);
        end
        for j = 1:size(k,2)
            K{j}(:,i - BurnIn)     = k(:,j);
            KA{j}(:,i - BurnIn)    = kA(:,j);
        end
    end
end

for j = 1:numel(beta) %To get the central tendency of the model. 
    Beta(:,j)  = prctile(beta{j}',50)';
    BetaA(:,j) = prctile(betaA{j}',50)';
end

E               = y - Beta(:,1);
K               = E'*Beta(:,2:end)*pinv(Beta(:,2:end)'*Beta(:,2:end));
K               = max(abs(prctile(K,[0.5 99.5]')));
Beta(:,2:end)   = Beta(:,2:end).*K; %To scale the coefficients so that K = +/- 1 represents extreme values. 

E               = yA - BetaA(:,1);
KA              = E'*BetaA(:,2:end)*pinv(BetaA(:,2:end)'*BetaA(:,2:end));
KA              = max(abs(prctile(KA,[0.5 99.5]')));
BetaA(:,2:end)  = BetaA(:,2:end).*KA;

for j = 2:numel(beta)
    beta{j}    = beta{j}*K(j - 1);
    betaA{j}   = betaA{j}*KA(j - 1);
end

hold off
coloR           = {[1 1 1],[0 0 1],[1 0 0],[0 1 0],[1 1 0]};
k               = 1;
s               = 200; %To define the number of samples for a plot. 
plot(x,exp(beta{1}(:,1:s)),'color',[coloR{1} 0.05])
hold on
for j = 2:5
    for h = 1:2
        plot(x,exp(beta{1}(:,1:s) + beta{j}(:,1:s)*k*(-1)^h),'color',[coloR{j} 0.05]) %To change one parameter at the time, alternating the sign of k.
    end
end
