%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate CF Transtion - Actual Prices  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
EXP_KAPPA = 1;

% Simulation Parameters
simT = 49; % number of periods to simulate forward

load('rSS_pf')
load('L1SS_pf');
load('LU1SS_pf');
load('PSS_pf');
load('pKSS_pf');
load('wageSS_pf');

load('steadyInitialRun_pf')
ssC = load('steadyInitialRun_pf');
ssC.PSS = PSS_pf;
ssC.VSS = VPath(:,end);

PInit = PSS_pf;
wInit = wageSS_pf(:,end);
pKInit = pKSS_pf(end);
rInit = rSS_pf;

L0 = L1SS_pf;
LU0 = LU1SS_pf;
L1 = L0;
LU1 = LU0;

%%% Changes Term
% 1. Productivity (sensitive to capital price assumptions...)
simDELZ = zeros(NINDS,simT);
simDELZ(:,1:T-1) = DELZF(:,2:end);

% 2. Trade Prices 
simDELPF = zeros(NTRADABLES,simT);
simDELPF(:,1:T-1) = reshape(DELPF_D(DELPF_D(:,1)>=2,3),NTRADABLES,T-1);

%%% Levels Terms
% 1. Capital Income
simCAPINC = zeros(1,simT);
simCAPINC(1:T-1) = CAPINC_D(3:end,2);
simCAPINC(T:end) = simCAPINC(T-1);

% 2. Foreign Demand Term
simEXPSHIFT = zeros(NTRADABLES,simT);
simEXPSHIFT(:,1:T-1) = reshape(EXPSHIFT_D(EXPSHIFT_D(:,1)>=2,3),NTRADABLES,T-1);
simEXPSHIFT(:,T:end) = repmat(simEXPSHIFT(:,T-1),1,simT-T+1);

% Shutting down changes in commodities prices... (inds 2 and 11)
simDELPF(2,:) = 0;
simDELPF(11,:) = 0;

% One time shock AT TIME = 10!
simDELPF(:,10) = sum(simDELPF,2);
simDELPF(:,1:9) = 0;
simDELPF(:,11:end) = 0;


% Fixing stuff
delZ = zeros(size(delZ));
KINC = ssC.KINC;
EXPSHIFT = EXP_KAPPA*ssC.EXPSHIFT; % Initialize. Must be Corrected!
P = PInit;

%% Loop for Counterfactual
% Optimization Parameters
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);

% Initializing Vectors
wVecT = zeros(nOccs,simT+1);
wVecT(:,1) = wInit;

rVecT = zeros(length(rInit),simT+1);
rVecT(:,1) = rInit;

pKVecT = zeros(1,simT+1);
pKVecT(1) = pKInit;

% Other stuff to keep track of.
LPath = zeros(length(L0),simT);
LUPath = zeros(length(L0),simT);
PPath = zeros(simT+1,1);
gdpPath = zeros(simT+1,1);
UPath = zeros(simT,1);

VPath = zeros(length(L1),simT);


gdpPath(1) = (ssC.WINC + ssC.KINC*pKInit)/ssC.PSS;
LPath(:,1) = L1;
LUPath(:,1) = LU1;
PPath(1) = PInit;
UPath(1) = sum(LU1/1000);
VPath(:,1) = ssC.VSS;


sprintf('*** STARTING ***\n***INITIAL VARIABLES***\nPrice Level: %2.4g\nReal GDP: %2.4g\nUnemployment Rate: %2.4g',PPath(1),gdpPath(1),UPath(1)*100)

for j = 1:simT
    sprintf('\nStarting Round %d\n',j)
   
    % Step 1: Update Pop Distribution
    L0 = L1;
    LU0 = LU1;
    
    tic
    % Step 2: Update the current value of relative prices
    r0 = rVecT(:,j)';
    w0 = wVecT(:,j);
    pK0 = pKVecT(:,j);
    
    delPF = simDELPF(:,j);
    EXPSHIFT = EXPSHIFT.*exp((1-sigma)*delPF);
    
    sprintf('Mean Change in Import Price: %4.2g\n', mean(delPF))
    
    % Step 3: Update current wage guess ,delPK,delPF,delZ
    x1 = lsqnonlin(@(x) excessDemandCapital(x(1:nOccs),x(end),w0,pK0,P,L0,LU0,n,un,b,params,techParams,r0,delPF,zeros(size(delZ)),KINC,EXPSHIFT,TMAT,TMAT_U,cMat,hMat,state),[wVecT(:,j);pKVecT(j)],[],[],options);
    wVecT(:,j+1) = x1(1:nOccs);
    pKVecT(:,j+1) = x1(end);
    pK1 = pKVecT(:,j+1);
    
    % Step 3b: Update prices...
    %sprintf('\nBleep\n')
    
    % Change in Endogenous Prices
    delW  = wVecT(:,j+1) - wVecT(:,j);
    delPK = log(pKVecT(:,j+1)) - log(pKVecT(:,j));
     
    % Given new wages and the initial price calculate new price
    MIDMATBM0 = [bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]];
    
    % Calculates changes in relative prices of tradables and levels of prices of non-tradables
    delPrices = (eye(NINDS) - MIDMATBM0)^(-1) * (-delZ + [BTT - eye(NTRADABLES); BNT]*delPF + BK*delPK + BL*delW);

    % Updates relative prices
    r1 = r0.*exp(delPrices(1:NTRADABLES))';
    rVecT(:,j+1) = r1';
    
    % For tradables do NOT want change in relative prices but change in CES bundle (see paper for details)
    adjTerm = 1/(1-sigma)*(log(1+r1.^(1-sigma)) - log(1+r0.^(1-sigma)));
    delP = delPrices;
    delP(1:NTRADABLES) = delPF + adjTerm';
    P = P*exp(dot(ALPHA(:,2),delP));
    PPath(j+1) = P;
    
    % Step 4: Calculate new distribution of workers
    [HSt,piPDF, piPDF_U, VPath(:,j+1)] = laborSupplyQF_n(wVecT(:,j+1)-log(P),L0,LU0,n,un,b,params,TMAT,TMAT_U,cMat,hMat,state);
    
    %sprintf('\nBloop\n')
    
    L1 = zeros(size(L0));
    LU1 = zeros(size(L0));

    for a = 0:(nAge-2)
        for s = 0:(nTypes-1)
            for t = 0:(nTen-1)
                for o = 1:nOccs

                    dex = nTypes*nTen*nOccs*a + nTen*nOccs*s + nOccs*t;
                    L = L0(dex + o);
                    LU = LU0(dex + o);

                    % Employed
                    for op = 1:nOccs
                        % Same job (more tenure)
                        if o == op
                           L1(dex + nTypes*nTen*nOccs + nOccs*(t<(nTen-1)) + o) = L1(dex + nTypes*nTen*nOccs + nOccs*(t<(nTen-1)) + o) + L*piPDF(dex+o,o)+LU*piPDF_U(dex+o,o);
                        % New job (reset tenure)
                        else
                           L1(dex + nTypes*nTen*nOccs - nOccs*t + op) = L1(dex + nTypes*nTen*nOccs - nOccs*t + op) + L*piPDF(dex+o,op)+LU*piPDF_U(dex+o,op);
                        end
                    end

                     % Newly unemployed
                     LU1(dex + nTypes*nTen*nOccs + o) = LU1(dex + nTypes*nTen*nOccs + o) + L*(1-sum(piPDF(dex+o,:),2));

                     % Continuing unemployed (lose tenure)
                     LU1(dex + nTypes*nTen*nOccs - nOccs*t + o) = LU1(dex + nTypes*nTen*nOccs - nOccs*t + o) + LU*(1-sum(piPDF_U(dex+o,:),2));
                end
            end
        end
    end

    % Repopulating with Younglings
    popFactor = ((sum(L0) + sum(LU0)) - (sum(L1) + sum(LU1)))/(sum(L0(state(:,1)==1) + LU0(state(:,1)==1)));

    L1(state(:,1)==1) = popFactor*L0(state(:,1)==1);
    LU1(state(:,1)==1) = popFactor*LU0(state(:,1)==1);
    UPath(j+1) = sum(LU1/1000);
    
    LPath(:,j+1) = L1;
    LUPath(:,j+1) = LU1;
    
    WINC = dot(HSt,exp(wVecT(:,j+1)));
    gdpPath(j+1) = (WINC + KINC*pK1)/P;
    %gdpPath(j+1) = (WINC/P + KINC*pK1);
    
    % Puts income in current units of income
    INCOME  = WINC+KINC*pKVecT(:,j+1);

    EXPTERM = EXPSHIFT;
    %EXPTERM = P.^sigma*EXPSHIFT;
    %EXPTERM = EXPSHIFT.*P;
    
    % Calculating Labor Demand Items ,delPK,delPF,delZ
    [HDt, Et] = laborDemandF_real(techParams,r1,INCOME,EXPTERM,wVecT(:,j+1));
    
    % Imports
    MIDMATAL2 = ALPHA(1:NTRADABLES,2).*(1-r1'.^(1-sigma)./(1+r1'.^(1-sigma)));
    MIDMATBM2 = bsxfun(@times,[BTT;BNT],1-r1.^(1-sigma)./(1+r1.^(1-sigma)))';
    IMPORTS = (MIDMATBM2*Et + MIDMATAL2*(INCOME));
    EXPORTS = EXPTERM.*r1'.^(1-sigma);
    
    tElapse = toc;
    % Step 5: Calculate error
    error = norm(exp(wVecT(:,j+1))-exp(wVecT(:,j)),2)/norm(exp(wVecT(:,j)));
    
    % Step 6: Display
    sprintf('That Loop Took %0.5g Seconds\nIteration: %d\nError: %0.5g\nPopulation Measure is %4.1g\nNX = %2.4g\nPrice Level: %2.4g\nReal GDP: %2.4g\nUnemployment Rate: %2.4g',tElapse,j,error,sum(L1)+sum(LU1),(sum(EXPORTS) - sum(IMPORTS)),P,gdpPath(j+1),UPath(j+1)*100)
end

timingPaths.simDELPF = simDELPF;
timingPaths.simDELZ = simDELZ;
timingPaths.wagePathMyopia = wVecT;
timingPaths.pKPathMyopia = pKVecT;
timingPaths.PPathMyopia = PPath;
timingPaths.VPathMyopia = VPath;

save('timingPaths','timingPaths')
