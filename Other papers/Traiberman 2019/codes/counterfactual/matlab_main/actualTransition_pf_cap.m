%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATING PERFECT FORESIGHT STEADY STATE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
EXP_KAPPA = 1;

% Loading Exogenous Variables and Myopic Guess
load('actualPaths')
load('steadyInitialRun_pf')

delPFpath = actualPaths.simDELPF;  
delZpath =  zeros(size(actualPaths.simDELZ));

EXPSHIFTInit = EXP_KAPPA*EXPSHIFT_D(EXPSHIFT_D(:,1)==1,3);

% Steady State for Initialization
load('rSS_pf')
load('L1SS_pf');
load('LU1SS_pf');
load('PSS_pf');
load('pKSS_pf');
load('wageSS_pf');

L1 = L1SS_pf;
LU1 = LU1SS_pf;
w0 = wageSS_pf(:,end);

PInit = PSS_pf;
wInit = wageSS_pf(:,end);
pKInit = pKSS_pf(end);
rInit = rSS_pf;

% First Guess: Myopia
wPath0 = actualPaths.wagePathMyopia;
pKPath0 = actualPaths.pKPathMyopia;

% Preallocate some empty cells
numT = size(wPath0,2);

wPath1 = zeros(size(wPath0));
pKPath1 = zeros(size(pKPath0));
rPath = zeros(length(rInit),numT);

exportPath= zeros(NTRADABLES,numT);
importPath= zeros(NTRADABLES,numT);


% Algo Parameters
options = optimset('Display','iter','MaxFunEvals',10000,'MaxIter',1000);
outer_error = 1;
outer_tol = .0005;
loop = 0;

MAXLOOPS = 50;
MINLOOPS = 5;
GS_LAMBDA = .875;

tStart_full = tic;

delZ = zeros(size(delZpath(:,1)));

diary off
diary(strcat('actualTransitionLog',date,'_',num2str(hour(datetime('now'))),'.txt'))
diary on

%%
while outer_error>outer_tol && loop<MAXLOOPS;
    tStart = tic;
    loop = loop+1;
    
    sprintf('*** BEGINNING LOOP %d ***',loop)
    
    % --INITIAL PERIOD--
    % Price levels
    P0 = PInit;
    P = P0;
    r0 = rInit;
    EXPSHIFT = EXPSHIFTInit;
    
    % First period wages and capital prices
    wPath1(:,1) = wInit;
    pKPath1(1) = pKInit;
    rPath(:,1) = rInit';

    % Initial Population Measure
    L0 = L1SS_pf;
    LU0 = LU1SS_pf;
    L1 = L0;
    LU1 = LU0;

    % Step 0: Solve for price changes implied by wage and pK path
    tic
    Ppath = getPriceChangesCap(wPath0,pKPath0,techParams,P0,r0,delPFpath,delZpath);
    toc
    
    tic
    % Step 1: From {w}_o,t solve for {V}_o,t
    [VPath, VPath_U] = getContinuationValues_n(bsxfun(@plus,wPath0,-log(Ppath)),n,un,params,state,TMAT,TMAT_U,cMat,hMat);

    % Shift by one...
    VPath = [VPath(:,2:end) VPath(:,end)];
    VPath_U = [VPath_U(:,2:end) VPath_U(:,end)];

    toc

    % --LOOP FOR REMAINING PERIODS--
    for j = 2:numT;
        
        % Step 0: Update prices - pay attention to indexing
        delPF = delPFpath(:,j-1);  
        delZ = delZpath(:,j-1);
        EXPSHIFT = EXPSHIFT.*exp((1-sigma)*delPF);
        
        sprintf('Mean Change in Import Price: %4.2g\n', mean(delPF))
        
        % Step 1: Update Pop Distribution
        L0 = L1;
        LU0 = LU1;

        tic
        % Step 2: Update the current value of relative prices
        r0 = rPath(:,j-1)';

        % Step 3: Update current wage guess [delPK,delPF,delZ]
        %wPath1(:,j) = lsqnonlin(@(x) relExcessDemand_pf(x,VPath(:,j-1),VPath_U(:,j-1),wPath1(:,j-1),P,L0,LU0,n,un,b,params,techParams,r0,delPK,delPF,delZ,KINC,EXPSHIFT,TA,TA_U,cMat,hMat,state),wPath0(:,j),[],[],options);
        x1  = lsqnonlin(@(x) excessDemandCapital_pf(x(1:nOccs),x(end),wPath1(:,j-1),pKPath1(:,j-1),VPath(:,j-1),VPath_U(:,j-1),P,L0,LU0,n,un,b,params,techParams,r0,delPF,delZ,KINC,EXPSHIFT,TAMAT,TAMAT_U,cMat,hMat,state),[wPath0(:,j);pKPath0(:,j)],[],[],options);
        
        wPath1(:,j) = x1(1:nOccs);
        pKPath1(j) = x1(end);
        pK1 = pKPath1(j);
        
        % Step 4: Update prices
        % Change in Wages
        delW = wPath1(:,j) - wPath1(:,j-1);
        delPK = log(pKPath1(:,j)) - log(pKPath1(:,j-1));

        % Given new wages and the initial price calculate new price
        MIDMATBM0 = [bsxfun(@times,[BTT;BNT],r0.^(1-sigma)./(1+r0.^(1-sigma))) [BTN;BNN]];

        % Calculates changes in relative prices of tradables and levels of prices of non-tradables
        delPrices = (eye(NINDS) - MIDMATBM0)^(-1) * (-delZ + [BTT - eye(NTRADABLES); BNT]*delPF + BK*delPK + BL*delW);

        % Updates relative prices
        r1 = r0.*exp(delPrices(1:NTRADABLES))';
        rPath(:,j) = r1';

        % Adjustment term for relative prices
        adjTerm = 1/(1-sigma)*(log(1+r1.^(1-sigma)) - log(1+r0.^(1-sigma)));

        % For tradables do NOT want change in relative prices but change in CES bundle 
        delP = delPrices;
        delP(1:NTRADABLES) = delPF + adjTerm';
        
        % Finally update price with C-D prefs
        P = P*exp(dot(ALPHA(:,2),delP));

        % Step 4: Calculate new distribution of workers
        [HSt,piPDF, piPDF_U] = laborSupplyF_pf_new(wPath1(:,j)-log(P),VPath(:,j-1),VPath_U(:,j-1),L0,LU0,n,un,b,params,TAMAT,TAMAT_U,cMat,hMat,state);

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

        WINC = dot(HSt,exp(wPath1(:,j)));

        % Puts income in current units of currency
        %INCOME  = WINC+KINC*P;
        INCOME  = WINC+KINC*pK1;
        
        EXPTERM = EXPSHIFT;
        %EXPTERM = P.^sigma*EXPSHIFT;
        %EXPTERM = EXPSHIFT.*P;
        
        % Calculating Labor Demand Items ,delPK,delPF,delZ
        [HDt, Et] = laborDemandF_real(techParams,r1,INCOME,EXPTERM,wPath1(:,j));

        % Imports
        MIDMATAL2 = ALPHA(1:NTRADABLES,2).*(1-r1'.^(1-sigma)./(1+r1'.^(1-sigma)));
        MIDMATBM2 = bsxfun(@times,[BTT;BNT],1-r1.^(1-sigma)./(1+r1.^(1-sigma)))';
        IMPORTS = (MIDMATBM2*Et + MIDMATAL2*(INCOME));
        EXPORTS = EXPTERM.*r1'.^(1-sigma);
        
        exportPath(:,j) = EXPORTS;
        importPath(:,j) = IMPORTS;

        
        % Step 5: Calculate error
        error = norm(exp(wPath1(:,j))-exp(wPath0(:,j)),2)/norm(exp(wPath0(:,j)),2);
        
        % Step 6: Display
        tElapse = toc;
        sprintf('That Loop Took %0.5g Seconds\nIteration: %d\nError: %0.5g\nPopulation Measure is %4.1g\nNX = %2.4g\nPrice Level: %2.4g',tElapse,j,error,sum(L1)+sum(LU1),(sum(EXPORTS) - sum(IMPORTS)),P)
    end
    
    
    % UPDATE STEP
    %outer_error = norm(reshape((exp(wPath1)-exp(wPath0))./(1+exp(wPath0)),nOccs*numT,1),Inf);
    outer_error = norm(reshape(([exp(wPath1);pKPath1]-[exp(wPath0);pKPath0])./(1+[exp(wPath0);pKPath0]),nOccs*numT+numT,1),Inf);
    wPath0 = GS_LAMBDA*wPath0 + (1-GS_LAMBDA)*wPath1;
    pKPath0 = GS_LAMBDA*pKPath0 + (1-GS_LAMBDA)*pKPath1;
    
    sprintf('\n\n\n******FULL LOOP %d COMPLETE!******\nLoop Time: %2.6g seconds\nError:     %0.4g units\n\n\n',loop,toc(tStart_full),outer_error)
    
end


PpathT_pf = Ppath;
PT_pf = P;
pKT_pf = pKPath0;
rT_pf = r1;
wageT_pf = wPath0;
L1T_pf = L1;
LU1T_pf = LU1;
VPathT_pf = VPath;
VPathUT_pf = VPath_U;

save('actualRun_pf')
save('rT_pf','rT_pf')
save('L1T_pf','L1T_pf');
save('LU1T_pf','LU1T_pf');
save('PT_pf','PT_pf');
save('pKT_pf','pKT_pf');
save('wageT_pf','wageT_pf');
save('PpathT_pf','PpathT_pf');
save('VPathT_pf','VPathT_pf');
save('VPathUT_pf','VPathUT_pf');


sprintf('\n\n\n******PROCEDURE COMPLETE!******\nLoop Time:   %2.6g seconds\nFinal Error: %0.4g units\n\n\n',toc(tStart_full),outer_error)

diary off