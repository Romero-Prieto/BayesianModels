function [Age,b,ASFRc,country,cohort] = FertilityModel

pATh        = "/Users/lshjr3/Documents/Demographic_Methods/Fertility/";
lISt        = {'asfrTR'};
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

asfrTR      = OuT{1};
Age         = (13.25:0.50:55)';
s           = (asfrTR.Cohort == asfrTR.Year - asfrTR.Age);
asfrTR.T(s) = char('L'*ones(sum(s),1));
s           = (asfrTR.Cohort == asfrTR.Year - asfrTR.Age - 1);
asfrTR.T(s) = char('U'*ones(sum(s),1));
country     = tabulate(asfrTR.Code);
country     = country(:,1);
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


sEL         = {'AUT';'BEL';'CAN';'CHE';'CHL';'DEUTNP';'DNK';'ESP';'FIN';'FRATNP';'GBRTENW';'HUN';'IRL';'ITA';'JPN';'KOR';'NLD';'NOR';'POL';'PRT';'SWE';'TWN';'USA'};
f           = cell2mat(ASFRc(ismember(country,sEL)));
sET         = (~isnan(f(74,:)) & ~isnan(f(5,:)));
f           = recode(f(:,sET),NaN,0);
[~,chu]     = bspline(1:numel(Age),25,3,3,0);
F           = FertilityFitting(Age,f);
F           = exp(chu*log(F));

sET         = Age >= 13 & Age < 55;
Y           = log(F(sET,:)./sum(F(sET,:)))';
X           = ones(size(Y,1),1);
YY          = Y'*Y;
XX          = X'*X;
XY          = X'*Y;
b           = XX\XY;
EE          = YY - b'*XY - XY'*b + (b'*XX*b + (b'*XX*b)')/2;
[U,~,~]     = svd(EE);
b           = [b' ones(size(b')) U];
