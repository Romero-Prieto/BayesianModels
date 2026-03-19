function [Age,beta,ASFRc,country,cohort] = FertilityModel

pATh          = "/Users/lshjr3/Documents/mATLaB/FeRTiLItY/";
lISt          = {'asfrTR'};
for i = 1:numel(lISt)
    opt                   = detectImportOptions(char(pATh + string(lISt{i}) + ".txt"));
    opt.Delimiter         = {'\t'  ' ' '+' '-'};
    opt.DataLines         = [4,inf];
    opt.VariableNamesLine = 3;
    for j = 1:numel(opt.VariableTypes)
        if isequal(opt.VariableTypes{j},'char')
            opt.VariableTypes{j} = 'categorical';
        end
    end
    OuT{i}                = readtable(char(pATh + string(lISt{i}) + ".txt"),opt);
    clear opt j
end

asfrTR        = OuT{1};
Age           = (13.25:0.50:55)';
s             = (asfrTR.Cohort == asfrTR.Year - asfrTR.Age);
asfrTR.T(s)   = char('L'*ones(sum(s),1));
s             = (asfrTR.Cohort == asfrTR.Year - asfrTR.Age - 1);
asfrTR.T(s)   = char('U'*ones(sum(s),1));
country       = tabulate(asfrTR.Code);
country       = country(:,1);
for i = 1:numel(country)
    cohort{i}  = (min(asfrTR.Cohort(asfrTR.Code == country(i))):max(asfrTR.Cohort(asfrTR.Code == country(i))))';
    ASFRc{i}   = NaN(numel(Age),numel(cohort{i}));
    for j = 1:numel(cohort{i})
        sET                               = asfrTR(asfrTR.Code == country(i) & asfrTR.Cohort == cohort{i}(j),{'Year','Age','ASFR','Cohort','T'});
        sET.Age(sET.T == 'L')             = sET.Age(sET.T == 'L') + .25; 
        sET.Age(sET.T == 'U')             = sET.Age(sET.T == 'U') + .75;
        sET                               = sortrows(sET,'Age');
        ASFRc{i}(ismember(Age,sET.Age),j) = sET.ASFR/2;
    end
end


sEL           = {'AUT';'BEL';'CAN';'CHE';'CHL';'DEUTNP';'DNK';'ESP';'FIN';'FRATNP';'GBRTENW';'HUN';'IRL';'ITA';'JPN';'KOR';'NLD';'NOR';'POL';'PRT';'SWE';'TWN';'USA'};
f             = cell2mat(ASFRc(ismember(country,sEL)));
sET           = (~isnan(f(74,:)) & ~isnan(f(5,:)));
f             = recode(f(:,sET),NaN,0);
[~,chu]       = bspline(1:numel(Age),25,3,3,0);
F             = FertilityFitting(Age,f);
F             = exp(chu*log(F));

Z             = 0;
sET           = Age >= 15 & Age < 50;
H             = log(sum(F(sET,:)));
X             = H.^((0:Z)');
y             = log(F) - H;
R             = rank(y*y');
[u,s,v]       = svds(y,R);
U             = u;
Y             = (v*s)';             

YY            = Y*Y';
XX            = X*X';
XY            = X*Y';
theta         = XX\XY;
EE            = YY - XY'*theta - theta'*XY + theta'*XX*theta;
[u,~,~]       = svds(U*EE*U',R);
Uols          = u;

[ages,tables] = size(Y);
I             = sparse(eye(ages));
priorB        = {zeros(1 + Z,ages),ones(1 + Z,ages)*10^-4};
priorPsi      = {ones(ages,1)*10^-2,ones(ages,1)*10^-2};
priorT        = {triu(NaN(ages)) + tril(zeros(ages),-1),triu(NaN(ages)) + tril(ones(ages),-1)*10^-2};
G             = sum(~isnan(priorT{1}));
G             = G(G > 0);
lisT          = find(~isnan(priorT{2}));

C             = sparse(zeros(ages^2,numel(lisT)));
for i = 1:numel(lisT)
    C(lisT(i),i) = 1;
end
B             = theta;
Psi           = sparse(diag(gamrnd(priorPsi{1},priorPsi{2})));
T             = eye(ages);

IterAT        = 5000;
BurnIn        = 250;
for i = 1:BurnIn + IterAT
    EE       = YY - XY'*B - B'*XY + B'*XX*B;
    for j = 1:ages
        Psi(j,j)   = gamrnd(priorPsi{1}(j) + (tables)/2,1/(priorPsi{2}(j) + 1/2*T(:,j)'*EE*T(:,j)));      
    end        
    
    L        = mat2cell(C'*kron(Psi,EE)*C + sparse(diag(priorT{2}(lisT))),G,G);
    V        = mat2cell(zeros(size(C,2)),G,G);
    for j = 1:numel(G)
        V{j,j} = sparse(pinv(full(L{j,j})));
        L{j,j} = chol(L{j,j})';        
    end

    L        = cell2mat(L);
    V        = cell2mat(V);
    Mu       = -C'*kron(Psi,EE)*reshape(I,[],1) + sparse(diag(priorT{2}(lisT))*priorT{1}(lisT));
    T(lisT)  = V*Mu + L'\normrnd(0,1,numel(lisT),1);
    Phi      = T*Psi*T';
    clear L V Mu

    V        = kron(Phi,XX) + diag(reshape(priorB{2},[],1));
    L        = chol(sparse(V))';
    V        = pinv(V);
    Mu       = kron(Phi,XX)*reshape(theta,[],1) + diag(reshape(priorB{2},[],1))*reshape(priorB{1},[],1);
    B        = V*Mu + L'\normrnd(0,1,numel(B),1);
    B        = reshape(B',Z + 1,ages);
    clear L V Mu

    [u,~,~]  = svds(U*pinv(full(Phi))*U',R);
    u        = u.*sign(sum(Uols.*u));
    u        = [U*B' u];
    k        = (y - u(:,1))'*u(:,2:end);

    if i > BurnIn
        for j = 1:size(u,2)
            beta{j}(:,i - BurnIn) = u(:,j);
        end
        K = [K;k];
    else
        K = [];
    end
end
