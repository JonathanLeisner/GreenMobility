%%%%%%%
% Calculate Steady State
%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithm:
% 1. Initialize from 1996 to 1997 using realized initial wages and realized
%    changes in productivity, capital prices, etc.
% 
% 2. Assume no further changes from 1996 to 1997 and iterate to steady
%    state on distribution of workers.
% 
% 3. In each period update new entrants from distribution of youngest
%    entrants in base year.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Optional Reset and Reload of Data
clear;
clc;

readInDemandData
readInParameters

%% PARAMETERIZATION%%
% Model Parameters
nTypes   =  6;
nTen     =  6;
nAge     = 35;
nOccs    =  38;
nChars   =  size(distMatrix,2);
sizeV    =  nTypes*nTen*nOccs*nAge;

% State = (age, type, tenure, occ)
state = INITIALEDF(:,2:5);

% Mobility Parameters
rho = RHO;
uCost = -GAMMA(2) * RHO;

fAge_o = SEARCH(1);
fAge2_o = SEARCH(2);
fType_o = SEARCH(3:end);

% Utility Parameters
eta  = ETA(:) * RHO;
beta = .96;

% Mincer Parameters
bAge_o  = MINCER(MINCER(:,1)==1 & MINCER(:,2)==1,4);
bAge2 = MINCER(MINCER(:,1)==1 & MINCER(:,2)==1,5);
bTen  = MINCER(MINCER(:,1)==1 & MINCER(:,2)==1,6);
bRMSE  = MINCER(MINCER(:,1)==1 & MINCER(:,2)==1,7);
bType = reshape(MINCER(MINCER(:,1)==1,end),nTypes,nOccs)';

bType(bType(:,1:nTypes-1)==0)=-10000;
bType(bType(:,5)==-10000,6) = -10000;

bType(bType(:,3)~=0 & bType(:,3)~=-10000 & bType(:,4)==-10000,4) = 0;

% Unemployment Parameters
uAge_o    = NONEMP(1);
uAge2_o   = NONEMP(2);
uCons_o   = NONEMP(3:end);


% Wages
w = MINCER(MINCER(:,1)==2 & MINCER(:,2)==1,end-1);


% Adjusting for age from [25-59] to [0,34]
% 1. Wage Parameters
%w       = w + 25*bAge_o + 625*bAge2 + bRMSE.^2/2;
w       = w + 25*bAge_o + 625*bAge2;
bAge    = bAge_o + 50*bAge2;

% 2. Search Parameters
%fType   = fType_o + 25*fAge_o + 625*fAge2_o;
%fAge    = fAge_o + 50*fAge2_o;
%fAge2   = fAge2_o;
% Corrections already made in estimation phase
fType = fType_o;
fAge = fAge_o;
fAge2 = fAge2_o;


% 3. Unemployment Parameters
%uCons   = uCons_o + 25*uAge_o + 625*uAge2_o;
%uAge    = uAge_o + 50*uAge2_o;
%uAge2   = uAge2_o;
uCons = uCons_o;
uAge = uAge_o;
uAge2 = uAge2_o;

%% LOAD STRUCTS FOR FUNCTIONS %%

%1. Model Parameters
n.Types	=	nTypes;
n.Ten	=	nTen;
n.Age	=	nAge;
n.Occs	=	nOccs;

% 2. Discounts/Variances
params.rho	=	rho;
params.beta	=	beta;

% 3. Stage Payout Parameters
u.eta	=	eta;

% 4. Unemployment Parameters
un.Cost	=	uCost;
un.Age	=	uAge;
un.Age2	=	uAge2;
un.Cons	=	uCons;

% 5. Search Function
f.age	=	fAge;
f.age2	=	fAge2;
f.type	=	fType;

% 6. Wage Function
b.age	=	bAge;
b.age2	=	bAge2;
b.Ten	=	bTen;
b.RMSE	=	bRMSE;
b.Type	=	bType;


% 7. Costs
costMat = zeros(nOccs, nOccs);
for i = 1:nOccs
    for j = 1:nOccs
        costMat(i,j) = -RHO*(exp(distMatrix(nOccs*(i-1)+j,3:end)*GAMMA(3:end) + GAMMA(1)))*(j~=i);
    end
end

costs = cell(1,nTypes);

for j = 1:nTypes
    costs{j} = costMat;
    costs{j}(:,bType(:,j)'==-10000)=-10000;
end

% 8. Labor Demand Parameters
techParams.NINDS = NINDS;
techParams.NTRADABLES = NTRADABLES;

% IO Matrix
techParams.BTT = BTT;
techParams.BTN = BTN;
techParams.BNT = BNT;
techParams.BNN = BNN;

techParams.BL = BL;
techParams.BK = BK;

% Consumer Parameters
techParams.ALPHA = ALPHA;
techParams.sigma = sigma;

%% Pre-Allocating Some Matrices %%
% State transitions
TAMAT = zeros(nTypes*nTen*nOccs*nAge,nOccs);
TAMAT_U = zeros(nTypes*nTen*nOccs*nAge,1);
for a = 0:(nAge-1)
    for s = 0:(nTypes-1)
        for t = 0:(nTen-1)
            for o = 1:nOccs
                %j = nTen*nOccs*s + nOccs*t + o;
                j = nTypes*nTen*nOccs*a + nTen*nOccs*s + nOccs*t + o;

                TAMAT_U(j) = nTypes*nTen*nOccs*a + nTen*nOccs*s + o;  % Reset to 0 tenure if unemployed for >1 period

                for op = 1:nOccs
                    if o == op
                        TAMAT(j,op) = nTypes*nTen*nOccs*a + nTen*nOccs*s + nOccs*min(t+1,nTen-1) + op;  % Tenure up if you stay
                    else
                        TAMAT(j,op) = nTypes*nTen*nOccs*a + nTen*nOccs*s + op;                          % Else tenure resets
                    end
                end
            end
        end
    end
end

TMAT = zeros(nTypes*nTen*nOccs,nOccs);
TMAT_U = zeros(nTypes*nTen*nOccs,1);
for s = 0:(nTypes-1)
    for t = 0:(nTen-1)
        for o = 1:nOccs
            j = nTen*nOccs*s + nOccs*t + o;
            
            TMAT_U(j) = nTen*nOccs*s + o;  % Reset to 0 tenure if unemployed for >1 period
            
            for op = 1:nOccs
                if o == op
                    TMAT(j,op) = nTen*nOccs*s + nOccs*min(t+1,nTen-1) + op;  % Tenure up if you stay
                else
                    TMAT(j,op) = nTen*nOccs*s + op;                          % Else tenure resets
                end
            end
        end
    end
end

% New Stage Payoff Method
cMat = zeros(sizeV,nOccs);
hMat = zeros(sizeV,nOccs);
j = 0;
for a = 0:(nAge-1)
    for s = 1:nTypes
        for t = 1:nTen
            for o = 1:nOccs
                j = j+1;
                % No benefit for youngs (initial conditions issues...)
                %cMat(j,:) = exp(fAge2*a^2 + fAge*a + fType(s))*costs{s}(o,:) + eta';
                
                % Cheap for youngs
                cMat(j,:) = (.1*(a==0) + 1*(a>0))*exp(fAge2*a^2 + fAge*a + fType(s))*costs{s}(o,:) + eta';
                
                % Sets non-existent pairs to super high costs
                cMat(j,costs{s}(o,:)<-5000) = -10000;
            end
        end
    end
end


j = 0;
A = eye(nOccs);
for a = 0:(nAge-1)
    for s = 1:nTypes
        for t = 1:nTen
            for o = 1:nOccs
                j = j+1;
                hMat(j,:) = exp(bAge*a + bAge2*a^2 + bType(:,s) + bTen.*A(:,o)*t + bRMSE.^2/2)';
            end
        end
    end
end

%% INITIALIZATION %%
% Initial Population Measure
L0 = INITIALEDF(:,6);
LU0 = INITIALEDF_U (:,6);
L1 = L0;
LU1 = LU0;

% Exogeneous Price and Productivity Shifts
period = 1;
delPK = DELR_D(DELR_D(:,1)==period,2);
delPF = DELPF_D(DELPF_D(:,1)==period,3);  
delZ = DELZF(:,period);
delDA = DELRELDEMAND_D(DELRELDEMAND_D(:,1)==period,3);
EXPSHIFT = EXPSHIFT_D(EXPSHIFT_D(:,1)==period,3);

rInit = RELPRICES_D(RELPRICES_D(:,1)==period-1,3)';

wInit = MINCER(MINCER(:,1)==period & MINCER(:,2)==1,end-1);
%wInit = wInit + 25*bAge_o + 625*bAge2 + bRMSE.^2/2;
wInit = wInit + 25*bAge_o + 625*bAge2;

% Exogeneous Capital Stock
KINC = CAPINC_D(CAPINC_D(:,1)==period,2);

% Initialize the price level
P = 1;

%% Loop for OTHER PERIODS
% Programming Parameters
maxIters = 50;
tol = .0025; % A low tolerance is needed because heavy history-dependence implies slow convergence
error = 1;
j = 1;
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);


% Initializing Vectors
wVec = zeros(nOccs,maxIters);
wVec(:,1) = wInit;

rVec = zeros(length(rInit),maxIters);
rVec(:,1) = rInit;

pKVec = zeros(1,maxIters);
pKVec(1) = 1;

VSS = zeros(sizeV,1);

% Stuff that's zero...
delPF = zeros(size(delPF));
delZ = zeros(size(delZ));

while error>tol && j<maxIters+1
    % Step 0: Update Iteration Count
    j = j + 1;
    
    % Step 1: Update Pop Distribution
    L0 = L1;
    LU0 = LU1;
    
    tic
    % Step 2: Update the current value of relative prices
    r0 = rVec(:,j-1)';
    w0 = wVec(:,j-1);
    pK0 = pKVec(:,j-1);
    
    % Step 3: Update current wage guess ,delPK,delPF,delZ
    %wVec(:,j) = lsqnonlin(@(x) relExcessDemand_real_n(x,wVec(:,j-1),P,L0,LU0,n,un,b,params,techParams,r0,zeros(size(delPK)),zeros(size(delPF)),zeros(size(delZ)),KINC,EXPSHIFT,TMAT,TMAT_U,cMat,hMat,state),wVec(:,j-1),[],[],options);
    x1 = lsqnonlin(@(x) excessDemandCapital(x(1:nOccs),x(end),w0,pK0,P,L0,LU0,n,un,b,params,techParams,r0,zeros(size(delPF)),zeros(size(delZ)),KINC,EXPSHIFT,TMAT,TMAT_U,cMat,hMat,state),[w0;pK0],[],[],options);
    wVec(:,j) = x1(1:nOccs);
    pKVec(:,j) = x1(end);
    pK1 = pKVec(:,j);
    
    % Step 3b: Update prices...
      % Change in Wages
    delW  = wVec(:,j) - wVec(:,j-1);
    delPK = log(pKVec(:,j)) - log(pKVec(:,j-1));
     
      % Given new wages and the initial price calculate new price
    MIDMATBM0 = [bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]];
    
      % Calculates changes in relative prices of tradables and levels of
      % prices of non-tradables
    delPrices = (eye(NINDS) - MIDMATBM0)^(-1) * (-delZ + [BTT - eye(NTRADABLES); BNT]*delPF + BK*delPK + BL*delW);

      % Updates relative prices
    r1 = r0.*exp(delPrices(1:NTRADABLES))';
    rVec(:,j) = r1';
    
      % Adjustment term for relative prices
    adjTerm = 1/(1-sigma)*(log(1+r1.^(1-sigma)) - log(1+r0.^(1-sigma)));
    delP = delPrices;
    
      % For tradables do NOT want change in relative prices but change in
      % CES bundle (see paper for details)
    delP(1:NTRADABLES) = delPF + adjTerm';
    P = P*exp(dot(ALPHA(:,2),delP));
    
    % Step 4: Calculate new distribution of workers
    [HSt,piPDF, piPDF_U, VSS] = laborSupplyQF_n(wVec(:,j)-log(P),L0,LU0,n,un,b,params,TMAT,TMAT_U,cMat,hMat,state);
    
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
    
    WINC = dot(HSt,exp(wVec(:,j)));
    
    % Puts income in current units of current
    INCOME  = WINC+KINC*pKVec(:,j);

    EXPTERM = EXPSHIFT;
    %EXPTERM = P.^sigma*EXPSHIFT;
    %EXPTERM = EXPSHIFT.*P;
    
    % Calculating Labor Demand Items ,delPK,delPF,delZ
    [HDt, Et] = laborDemandF_real(techParams,r1,INCOME,EXPTERM,wVec(:,j));
    
    % Imports
    MIDMATAL2 = ALPHA(1:NTRADABLES,2).*(1-r1'.^(1-sigma)./(1+r1'.^(1-sigma)));
    MIDMATBM2 = bsxfun(@times,[BTT;BNT],1-r1.^(1-sigma)./(1+r1.^(1-sigma)))';
    IMPORTS = (MIDMATBM2*Et + MIDMATAL2*(INCOME));
    EXPORTS = EXPTERM.*r1'.^(1-sigma);

    
    tElapse = toc;
    % Step 5: Calculate error
    error = norm(exp(wVec(:,j))-exp(wVec(:,j-1)),2)/norm(exp(wVec(:,j-1)));
    
    % Step 6: Display
    sprintf('That Loop Took %0.5g Seconds\nIteration: %d\nError: %0.5g\nPopulation Measure is %4.1g\nNX = %2.4g\nPrice Level: %2.4g',tElapse,j,error,sum(L1)+sum(LU1),sum(EXPORTS) - sum(IMPORTS),P)
    
    if j == 2
        L97 = L1;
        LU97 = LU1;
    end


end

finalIter = j;


wSS = wVec(:,finalIter);
pKSS = pKVec(:,finalIter);
rSS = rVec(:,finalIter);
PSS = P;

wVecSS = wVec(:,1:finalIter);
pKVecSS = pKVec(:,1:finalIter);

LSS = L1;
LUSS = LU1;

piSS = piPDF;

save('steadyStateCapital')
save('wVecSS','wVecSS')
save('pKVecSS','pKVecSS')
save('EXPORTS_SS','EXPORTS')
save('IMPORTS_SS','IMPORTS')