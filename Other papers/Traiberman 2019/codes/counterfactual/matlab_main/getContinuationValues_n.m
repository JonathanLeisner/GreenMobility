function [VPath, VPath_U] = getContinuationValues_n(wPath,n,un,params,state,T,T_U,cMat,hMat)

%%%%%%%%%%
% Inputs:
% wPath = path of wages (real or nominal doesn't matter)
% n = model parameters
% u = utility paramters
% un = unemployment parameters
% f = mobility parameters
% b = wage parameters
% 

%% Prelminary Constants
% Mathematical Constants
PSI = -psi(1);

% Model Parameters
nTypes   =  n.Types;
nTen     =  n.Ten;
nAge     =  n.Age;
nOccs    =  n.Occs;
sizeV    =  nTypes*nTen*nOccs*nAge;

numT = size(wPath,2);

%%

% Mobility Parameters
rho = params.rho;

% Utility Parameters
beta = params.beta;

% Unemployment Parameters
uCost   = un.Cost;
uAge    = un.Age;
uAge2   = un.Age2;
uCons   = un.Cons;

VPath = zeros(sizeV,numT);
VPath_U = zeros(sizeV,numT);



%% Step 1: Solve the TERMINAL period by assuming steady state wages forever

% Set wage.
wG = wPath(:,end);

% Solve Worker's Problem
% 0. Correct nonemployment
uCons = uCons+beta*PSI;

% 1. Stage Payoff
stagePay = (cMat + repmat(exp(wG'),sizeV,1).*hMat)/rho+beta*PSI;

EMAX = zeros(sizeV,1);
EMAX_U = zeros(sizeV,1);

% 2. Terminal Age
D = exp(stagePay(sizeV - nTypes*nTen*nOccs + 1:end,:)- beta*PSI);
EMAX(sizeV - nTypes*nTen*nOccs + 1:end,1) = log(exp(uCons(state(sizeV - nTypes*nTen*nOccs + 1:end,2)) + uAge*nAge + uAge2*nAge^2) + sum(D,2));

D = exp(uCost + stagePay(sizeV - nTypes*nTen*nOccs + 1:end,:)- beta*PSI);
EMAX_U(sizeV - nTypes*nTen*nOccs + 1:end,1) = log(exp(uCons(state(sizeV - nTypes*nTen*nOccs + 1:end,2)) + uAge*nAge + uAge2*nAge^2) + sum(D,2));


% 3. Preceding Ages
for a = 2:nAge
    beg = sizeV - a*nTypes*nTen*nOccs + 1;
    fin = sizeV - (a-1)*nTypes*nTen*nOccs;
    
    E = EMAX(beg+nTypes*nTen*nOccs:fin+nTypes*nTen*nOccs);        % Continuation values for employment
    E_U = EMAX_U(beg+nTypes*nTen*nOccs:fin+nTypes*nTen*nOccs);    % Continuation values for unemployment
    D = exp(stagePay(beg:fin,:) + beta*E(T));                     % Numerators
    
    EMAX(beg:fin) = log(exp(uCons(state(beg:fin,2)) + uAge*(nAge-a+1) + uAge2*(nAge-a+1)^2 + beta*E_U) + sum(D,2));         % New Continuation Values
    
    D = exp(uCost + stagePay(beg:fin,:) + beta*E(T));                        % Numerators if unemployed
    
    EMAX_U(beg:fin) = log(exp(uCons(state(beg:fin,2)) + uAge*(nAge-a+1) + uAge2*(nAge-a+1)^2 + beta*E_U(T_U)) + sum(D,2));  % New Unemp. Continuation Values
end

VPath(:,end) = EMAX;
VPath_U(:,end) = EMAX_U;

%% Now go backwards to solve remaining value functions...

% Set wage
%period = numT-1;

for period = numT-1:-1:1
    wG = wPath(:,period);

    % 1.  Update Stage Pay
    stagePay = (cMat + repmat(exp(wG'),sizeV,1).*hMat)/rho+beta*PSI;

    EMAX = zeros(sizeV,1);
    EMAX_U = zeros(sizeV,1);

    % 2. Terminal Age
    D = exp(stagePay(sizeV - nTypes*nTen*nOccs + 1:end,:)- beta*PSI);
    EMAX(sizeV - nTypes*nTen*nOccs + 1:end,1) = log(exp(uCons(state(sizeV - nTypes*nTen*nOccs + 1:end,2)) + uAge*nAge + uAge2*nAge^2) + sum(D,2));

    D = exp(uCost + stagePay(sizeV - nTypes*nTen*nOccs + 1:end,:)- beta*PSI);
    EMAX_U(sizeV - nTypes*nTen*nOccs + 1:end,1) = log(exp(uCons(state(sizeV - nTypes*nTen*nOccs + 1:end,2)) + uAge*nAge + uAge2*nAge^2) + sum(D,2));


    % 3. Preceding Ages
    for a = 2:nAge
        beg = sizeV - a*nTypes*nTen*nOccs + 1;
        fin = sizeV - (a-1)*nTypes*nTen*nOccs;

        E = VPath(beg+nTypes*nTen*nOccs:fin+nTypes*nTen*nOccs,period+1);        % Continuation values for employment using TOMORROW'S VALUES
        E_U = VPath_U(beg+nTypes*nTen*nOccs:fin+nTypes*nTen*nOccs,period+1);    % Continuation values for unemployment using TOMORROW'S VALUES
        D = exp(stagePay(beg:fin,:) + beta*E(T));                               % Numerators

        EMAX(beg:fin) = log(exp(uCons(state(beg:fin,2)) + uAge*(nAge-a+1) + uAge2*(nAge-a+1)^2 + beta*E_U) + sum(D,2));         % New Continuation Values

        D = exp(uCost + stagePay(beg:fin,:) + beta*E(T));                        % Numerators if unemployed

        EMAX_U(beg:fin) = log(exp(uCons(state(beg:fin,2)) + uAge*(nAge-a+1) + uAge2*(nAge-a+1)^2 + beta*E_U(T_U)) + sum(D,2));  % New Unemp. Continuation Values
    end

    VPath(:,period) = EMAX;
    VPath_U(:,period) = EMAX_U;
end