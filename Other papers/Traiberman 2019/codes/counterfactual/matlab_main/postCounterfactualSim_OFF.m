%figurePath = '/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/figures/';
%outputPath = '/Users/sharontraiberman/Dropbox/DanishLEEDProject/structuralEstimation_Final/counterfactualOutputs/';

figurePath = '~/Desktop/';
outputPath = '~/Desktop/';

load('steadyInitialRun_pf');
cf_pf = load('counterfactualRun_pf');
off_pf = load('OFF_counterfactualRun_pf');

newcode2disco;
L1SS_pf = off_pf.L1SS_pf;
LU1SS_pf = off_pf.L1SS_pf;
state = off_pf.state;
numT = off_pf.numT;

%%
%%%%%%%%%%%%%%%%%%%%%
% Preliminary Stuff %
%%%%%%%%%%%%%%%%%%%%%
close all
capT = 12; % number of periods to plot out

realWageT = bsxfun(@plus,off_pf.wagecf_pf',-log(off_pf.Ppathcf_pf'));
realWageCF = bsxfun(@plus,cf_pf.wagecf_pf',-log(cf_pf.Ppathcf_pf'));

tic
[LPathT, LPathUT, occDistPathT, piPathT, piUPathT] = makeLPath(L1SS_pf,LU1SS_pf,realWageT',off_pf.VPathcf_pf,off_pf.VPathUcf_pf,state,n,b,un,params,TAMAT,TAMAT_U,cMat,hMat,off_pf.numT);
toc

tic
[LPathCF, LPathUCF, occDistPathCF, piPathCF, piUPathCF] = makeLPath(L1SS_pf,LU1SS_pf,realWageCF',cf_pf.VPathcf_pf,cf_pf.VPathUcf_pf,state,n,b,un,params,TAMAT,TAMAT_U,cMat,hMat,cf_pf.numT);
toc



%% 
%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Simulations   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

tic
[historyT, switchesT, totalsT] = popSimulator(cumsum(L1SS_pf)/sum(L1SS_pf), piPathT,piUPathT, state, 100000, 479);
[historyCF, switchesCF, totalsCF] = popSimulator(cumsum(L1SS_pf)/sum(L1SS_pf), piPathCF,piUPathCF, state, 100000, 479);
toc

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Simulation Analysis   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From cells to matrices
historyTMat = cell2mat(historyT);
historyTMat(historyTMat(:,1)==0,:)=[];

historyCFMat = cell2mat(historyCF);
historyCFMat(historyCFMat(:,1)==0,:)=[];

wageTMat = zeros(length(historyTMat),1);
wageCFMat = zeros(length(historyCFMat),1);

tic
for j = 1:length(historyTMat)
    wageTMat(j) = hMat(historyTMat(j,3),historyTMat(j,7))*exp(realWageT(historyTMat(j,2),historyTMat(j,7))).*historyTMat(j,8);
    wageCFMat(j) = hMat(historyCFMat(j,3),historyCFMat(j,7))*exp(realWageCF(historyCFMat(j,2),historyCFMat(j,7))).*historyCFMat(j,8);
end
toc



% Writing to csv file for Stata
dlmwrite(strcat(outputPath,'historyT_OFF.csv'),[historyTMat wageTMat]);
dlmwrite(strcat(outputPath,'historyCF_OFF.csv'),[historyCFMat wageCFMat]);


%%
%%%%%%%%%%%%%%%%%%%%
%%%   Figures   %%%%
%%%%%%%%%%%%%%%%%%%%
%%%
% 1 -13 = Man
% 14 - 25 = Serv
% 26 - 32 = FIRE
% 33 - 38 = Health
%%%

meanManChange = mean((off_pf.wagecf_pf(1:13,1:capT) - cf_pf.wagecf_pf(1:13,1:capT))*100);
meanServChange = mean((off_pf.wagecf_pf(14:25,1:capT) - cf_pf.wagecf_pf(14:25,1:capT))*100);
meanFIREChange = mean((off_pf.wagecf_pf(26:32,1:capT) - cf_pf.wagecf_pf(26:32,1:capT))*100);
meanHEChange = mean((off_pf.wagecf_pf(33:38,1:capT) - cf_pf.wagecf_pf(33:38,1:capT))*100);


fig1_rev = figure;
hold on
nomManChange = plot(-1:capT-2,meanManChange,'-','Color','blue','LineWidth',1.3);
nomServChange = plot(-1:capT-2,meanServChange,'--','Color',[1 1 1]*.5,'LineWidth',1.3);
nomFIREChange = plot(-1:capT-2,meanFIREChange,'-*','Color',[1 1 1]*.5,'LineWidth',1.3);
nomHEChange = plot(-1:capT-2,meanHEChange,'-.','Color',[1 1 1]*.5,'LineWidth',1.3);

ylim([-14,4])

legend('Manufacturing','Other Services','FIRE','Health & Educ.')
ylabel('% Change in Skill Prices')
xlabel('Periods After Shock')

hgexport(fig1_rev,strcat(figurePath,'nominalSectorWageDelta_OFF'))

% 1. Changes in Skill Prices Over Time [Should we do a table?]
% 1a. Nominal
fig1 = figure;
hold on
nomServLines = plot(repmat((-1:capT-2)',1,25),(off_pf.wagecf_pf(14:end,1:capT)' - cf_pf.wagecf_pf(14:end,1:capT)')*100,'-','Color',[1 1 1]*.85,'LineWidth',.9);
nomManLines = plot(repmat((-1:capT-2)',1,13),(off_pf.wagecf_pf(1:13,1:capT)' - cf_pf.wagecf_pf(1:13,1:capT)')*100,'-','Color','blue','LineWidth',1.3);

nomServGroup = hggroup;
nomManGroup = hggroup;


set(nomServLines,'Parent',nomServGroup)
set(nomManLines,'Parent',nomManGroup)


set(get(get(nomServGroup,'Annotation'),'LegendInformation'), 'IconDisplayStyle','on'); 
set(get(get(nomManGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 

legend('Non-Manufacturing','Manufacturing')
ylabel('% Change in Skill Prices')
xlabel('Periods After Shock')




%title('Nominal Skill Prices')
hold off
hgexport(fig1,strcat(figurePath,'nominalWagesDelta_OFF'))


% 1b. Real
fig2 = figure;
hold on
realServLines = plot(repmat((-1:capT-2)',1,25),(realWageT(1:capT,14:end) - realWageCF(1:capT,14:end))*100,'-','Color',[1 1 1]*.85,'LineWidth',.9);
realManLines = plot(repmat((-1:capT-2)',1,13),(realWageT(1:capT,1:13) - realWageCF(1:capT,1:13))*100,'-','Color','blue','LineWidth',1.3);

realServGroup = hggroup;
realManGroup = hggroup;

set(realServLines,'Parent',realServGroup)
set(realManLines,'Parent',realManGroup)

set(get(get(realServGroup,'Annotation'),'LegendInformation'), 'IconDisplayStyle','on'); 
set(get(get(realManGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); 
legend('Non-Manufacturing','Manufacturing')




ylabel('% Change in Skill Prices')
xlabel('Periods After Shock')
%title('Real Skill Prices')
hold off
hgexport(fig2,strcat(figurePath,'realWagesDelta_OFF'))





%% Real Income

gdpPathT_pf = zeros(1,numT);
gdpPathcf_pf = zeros(1,numT);

for t = 1:numT
    gdpPathT_pf(t) = (dot(occDistPathT(:,t),exp(off_pf.wagecf_pf(:,t))) + off_pf.pKcf_pf(t)*KINC)/off_pf.Ppathcf_pf(t);
    gdpPathcf_pf(t) = (dot(occDistPathCF(:,t),exp(cf_pf.wagecf_pf(:,t))) + cf_pf.pKcf_pf(t)*KINC)/cf_pf.Ppathcf_pf(t);
end

fig3 = figure
plot((log(gdpPathT_pf(2:end))-log(gdpPathcf_pf(2:end)))*100,'-k','LineWidth',2)
title('Real GDP Relative to Steady State Path')
xlabel('Periods')
ylabel('% Difference in GDP')
hgexport(fig3,strcat(figurePath,'GDPpath_OFF'))





%% Variance Decomposition [w/ trade equilibrium as weighting]
wagecf_pf = cf_pf.wagecf_pf;
Ppathcf_pf = cf_pf.Ppathcf_pf;

wageT_pf = off_pf.wagecf_pf;
PpathT_pf = off_pf.Ppathcf_pf;

VDT = zeros(2,numT);
VDcf = zeros(2,numT);

for t = 1:numT
    wTrade = exp(wageT_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(PpathT_pf(t)));
    wNoTrade = exp(wagecf_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(Ppathcf_pf(t)));
    
    sTrade = LPathT(:,t);
    sTrade = sTrade/sum(sTrade); % Normalizing to weights
    
    sNoTrade = LPathCF(:,t);
    sNoTrade = sNoTrade/sum(sNoTrade); % Normalizing to weights
        
    mu = dot(wTrade,sTrade);
    var = dot((wTrade - mu).^2,sTrade);
    muG = accumarray(A(:,1),wTrade.*sTrade)./accumarray(A(:,1),sTrade);
    varG = dot((wTrade - muG(A(:,1))).^2,sTrade);
    
    VDT(1,t) = var;
    VDT(2,t) = var - varG;
    
    mu = dot(wNoTrade,sNoTrade);
    var = dot((wNoTrade - mu).^2,sNoTrade);
    muG = accumarray(A(:,1),wNoTrade.*sNoTrade)./accumarray(A(:,1),sNoTrade);
    varG = dot((wNoTrade - muG(A(:,1))).^2,sNoTrade);
    
    VDcf(1,t) = var;
    VDcf(2,t) = var - varG;

end


VD2 = zeros(2,numT);
VLEVEL2 = zeros(2,numT);

for t = 1:numT
    wTrade = exp(wageT_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(PpathT_pf(t)));
    wNoTrade = exp(wagecf_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(Ppathcf_pf(t)));
    
    diff = wTrade - wNoTrade;
    
    sTrade = LPathT(:,t);
    sTrade = sTrade/sum(sTrade); % Normalizing to weights
    
        
    mu = dot(diff,sTrade);
    var = dot((diff - mu).^2,sTrade);
    muG = accumarray(A(:,1),diff.*sTrade)./accumarray(A(:,1),sTrade);
    varG = dot((diff - muG(A(:,1))).^2,sTrade);
    
    VD2(1,t) = var;
    VD2(2,t) = var -varG;
    
    mu = dot(wTrade,sTrade);
    var = dot((wTrade - mu).^2,sTrade);
    muG = accumarray(A(:,1),wTrade.*sTrade)./accumarray(A(:,1),sTrade);
    varG = dot((wTrade - muG(A(:,1))).^2,sTrade);
    
    VLEVEL2(1,t) = var;
    VLEVEL2(2,t) = var - varG;
end



VDocc = zeros(2,numT);
VLEVELocc = zeros(2,numT);

for t = 1:numT
    wTrade = exp(wageT_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(PpathT_pf(t)));
    wNoTrade = exp(wagecf_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(Ppathcf_pf(t)));
    
    diff = wTrade - wNoTrade;
    
    sTrade = LPathT(:,t);
    sTrade = sTrade/sum(sTrade); % Normalizing to weights
    
        
    mu = dot(diff,sTrade);
    var = dot((diff - mu).^2,sTrade);
    muG = accumarray(A(:,2),diff.*sTrade)./accumarray(A(:,2),sTrade);
    varG = dot((diff - muG(A(:,2))).^2,sTrade);
    
    VDocc(1,t) = var;
    VDocc(2,t) = var -varG;
    
    mu = dot(wTrade,sTrade);
    var = dot((wTrade - mu).^2,sTrade);
    muG = accumarray(A(:,2),wTrade.*sTrade)./accumarray(A(:,2),sTrade);
    varG = dot((wTrade - muG(A(:,2))).^2,sTrade);
    
    VLEVELocc(1,t) = var;
    VLEVELocc(2,t) = var - varG;
end



VDsector = zeros(2,numT);
VLEVELsector = zeros(2,numT);

for t = 1:numT
    % Levels Version
    wTrade =exp(wageT_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(PpathT_pf(t)));
    wNoTrade =exp(wagecf_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(Ppathcf_pf(t)));
    
    % Log Version
    %wTrade = (wageT_pf(state(:,4),t));
    %wNoTrade = (wagecf_pf(state(:,4),t));
    
    
    diff = wTrade - wNoTrade;
    
    sTrade = LPathT(:,t);
    sTrade = sTrade/sum(sTrade); % Normalizing to weights
    
        
    mu = dot(diff,sTrade);
    var = dot((diff - mu).^2,sTrade);
    muG = accumarray(A(:,3),diff.*sTrade)./accumarray(A(:,3),sTrade);
    varG = dot((diff - muG(A(:,3))).^2,sTrade);
    
    VDsector(1,t) = var;
    VDsector(2,t) = var-varG;
    
    mu = dot(wTrade,sTrade);
    var = dot((wTrade - mu).^2,sTrade);
    muG = accumarray(A(:,3),wTrade.*sTrade)./accumarray(A(:,3),sTrade);
    varG = dot((wTrade - muG(A(:,3))).^2,sTrade);
    
    VLEVELsector(1,t) = var;
    VLEVELsector(2,t) = var-varG;
    
end

s = 1;

lastPeriod=20;

fig5 = figure;
hold on
plot((s:lastPeriod),VD2(2,s:lastPeriod)./VD2(1,s:lastPeriod),'-k','LineWidth',2)
plot((s:lastPeriod),VDocc(2,s:lastPeriod)./VDocc(1,s:lastPeriod),'--k','LineWidth',2)
plot((s:lastPeriod),VDsector(2,s:lastPeriod)./VDsector(1,s:lastPeriod),':k','LineWidth',2)
hold off
h_xlab = xlabel('Year');
h_ylab =ylabel('Fraction of Variation Explained');
%h_title =title('Comparing Variance in Outcomes Across Equilibrium');
h_legend=legend('Occ x Sector','Occupations Only','Sectors Only');
set(h_legend,'FontSize',14);
%set(h_title,'FontSize',14);
set(h_xlab,'FontSize',12);
set(h_ylab,'FontSize',12);
hgexport(fig5,strcat(figurePath,'dynamicVarianceDecomposition_OFF'))








figure
hold on
plot((s:lastPeriod),VLEVEL2(2,s:lastPeriod)./VLEVEL2(1,s:lastPeriod),'-k','LineWidth',2)
plot((s:lastPeriod),VLEVELocc(2,s:lastPeriod)./VLEVELocc(1,s:lastPeriod),'--k','LineWidth',2)
plot((s:lastPeriod),VLEVELsector(2,s:lastPeriod)./VLEVELsector(1,s:lastPeriod),':k','LineWidth',2)
hold off
h_xlab = xlabel('Year');
h_ylab =ylabel('Fraction of Variation Explained');
h_title =title('Levels');
h_legend=legend('Occ x Sector','Occupations Only','Sectors Only');
set(h_legend,'FontSize',14);
set(h_title,'FontSize',14);
set(h_xlab,'FontSize',12);
set(h_ylab,'FontSize',12);


levelMean_os = mean(VLEVEL2(2)./VLEVEL2(1));
levelMean_o = mean(VLEVELocc(2)./VLEVELocc(1));
levelMean_s =  mean(VLEVELsector(2)./VLEVELsector(1));

diffSRmean_os = VD2(2,2)./VD2(1,2);
diffSRmean_o = VDocc(2,2)./VDocc(1,2);
diffSRmean_s = VDsector(2,2)./VDsector(1,2);

diffLRmean_os = VD2(2,end)./VD2(1,end);
diffLRmean_o = VDocc(2,end)./VDocc(1,end);
diffLRmean_s = VDsector(2,end)./VDsector(1,end);

varDecompM_OFF = [levelMean_os levelMean_o levelMean_s;
diffSRmean_os diffSRmean_o diffSRmean_s;
diffLRmean_os diffLRmean_o diffLRmean_s]'*100;


dlmwrite(strcat(outputPath,'varDecomp_OFF.csv'),varDecompM_OFF);
%% - Skill Prices only (not income)
spVariance = zeros(1,numT);
spVariance_occ = zeros(1,numT);
spVariance_sec = zeros(1,numT);

diffVariance = zeros(1,numT);
diffVariance_occ = zeros(1,numT);
diffVariance_sec = zeros(1,numT);

AA = A(1:38,:);

% Skill Price Variance Decomposition
for t = 1:numT
    % Step 1: Weights are from the trade equilibrium
    weights = accumarray(state(:,4),LPathT(:,t))/sum(accumarray(state(:,4),LPathT(:,t)));
    
    incTrade = exp(wageT_pf(:,t))/PpathT_pf(t);
    incCF = exp(wagecf_pf(:,t))/Ppathcf_pf(t);
    diff = incTrade - incCF;
    
    % Step 2: Means for levels
    muT = dot(incTrade,weights);
    
    muGT_occ = accumarray(AA(:,2),weights.*incTrade)./accumarray(AA(:,2),weights);
    muGT_sec = accumarray(AA(:,3),weights.*incTrade)./accumarray(AA(:,3),weights);

    % Step 3: Variances for levels
    varT = dot((incTrade - mu).^2,weights);

    varGT_occ = dot((incTrade - muGT_occ(AA(:,2))).^2,weights);
    varGT_sec = dot((incTrade - muGT_sec(AA(:,3))).^2,weights);
    
    
    % Step 4: Means for diffs
    muD = dot(diff,weights);
    
    muGD_occ = accumarray(AA(:,2),weights.*diff)./accumarray(AA(:,2),weights);
    muGD_sec = accumarray(AA(:,3),weights.*diff)./accumarray(AA(:,3),weights);

    % Step 5: Variances for diffs
    varD = dot((diff - muD).^2,weights);

    varGD_occ = dot((diff - muGD_occ(AA(:,2))).^2,weights);
    varGD_sec = dot((diff - muGD_sec(AA(:,3))).^2,weights);

    % Step 6: Allocate to vectors
    spVariance(t)=varT;
    spVariance_occ(t)=varT - varGT_occ;
    spVariance_sec(t)=varT - varGT_sec;
    
    diffVariance(t)=varD;
    diffVariance_occ(t)=varD - varGD_occ;
    diffVariance_sec(t)=varD - varGD_sec;

end

figure
hold on
plot((s:lastPeriod),spVariance_occ(1,s:lastPeriod)./spVariance(1,s:lastPeriod),'--k','LineWidth',2)
plot((s:lastPeriod),spVariance_sec(1,s:lastPeriod)./spVariance(1,s:lastPeriod),':k','LineWidth',2)
hold off
h_xlab = xlabel('Year');
h_ylab =ylabel('Fraction of Variation Explained');
h_title =title('Levels');
h_legend=legend('Occupations Only','Sectors Only');
set(h_legend,'FontSize',14);
set(h_title,'FontSize',14);
set(h_xlab,'FontSize',12);
set(h_ylab,'FontSize',12);


figure
hold on
plot((s:lastPeriod),diffVariance_occ(1,s:lastPeriod)./diffVariance(1,s:lastPeriod),'--k','LineWidth',2)
plot((s:lastPeriod),diffVariance_sec(1,s:lastPeriod)./diffVariance(1,s:lastPeriod),':k','LineWidth',2)
hold off
h_xlab = xlabel('Year');
h_ylab =ylabel('Fraction of Variation Explained');
h_title =title('Differences');
h_legend=legend('Occupations Only','Sectors Only');
set(h_legend,'FontSize',14);
set(h_title,'FontSize',14);
set(h_xlab,'FontSize',12);
set(h_ylab,'FontSize',12);


%% 
t = 2;
wTrade = exp(wageT_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(PpathT_pf(t)));
wNoTrade = exp(wagecf_pf(state(:,4),t) + bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2))) - log(Ppathcf_pf(t)));
sTrade = LPathT(:,t);
sTrade = sTrade/sum(sTrade); % Normalizing to weights
mu_occ = accumarray(A(:,2),(wTrade-wNoTrade).*sTrade)./accumarray(A(:,2),sTrade);
mu_sec = accumarray(A(:,3),(wTrade-wNoTrade).*sTrade)./accumarray(A(:,3),sTrade);
mu_nohet = accumarray(A(:,1),(wTrade-wNoTrade).*sTrade)./accumarray(A(:,1),sTrade);

std(wTrade-wNoTrade)
std(mu_sec)
std(mu_occ)