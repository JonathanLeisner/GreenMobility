%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to solve the labor supply of workers given
% initial conditions and an assumption of myopia.
% 
% With myopic preferences the worker's problem is solved
% period by period. With EMAX equation,
% 
% EMAX_o = gamma + log(sum(exp(u_oo' + beta * EMAX_o')))
% 
% Worker state in this case is AGE, TENURE, OCCUPATION, TYPE
% 
% Indexing:
% Age  = {0,...,nAge-1}
% Ten  = {0,...,nTen-1}
% Type = {0,...,nType-1}
% Occ  = {1,...,nOccs}
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [HS, piPDF, piPDF_U, EMAX] = laborSupplyQF_n(w,L0,LU0,n,un,b,params,T,T_U,cMat,hMat,state)

%% Prelminary Constants
% Mathematical Constants
PSI = -psi(1);

% Model Parameters
nTypes   =  n.Types;
nTen     =  n.Ten;
nAge     =  n.Age;
nOccs    =  n.Occs;
sizeV    =  nTypes*nTen*nOccs*nAge;

%%

% Mobility Parameters
rho = params.rho;

% Utility Parameters
beta = params.beta;

% Mincer Parameters
bAge  = b.age;
bAge2 = b.age2;
bTen  = b.Ten;
bType = b.Type;
bRMSE = b.RMSE;

% Unemployment Parameters
uCost   = un.Cost;
uAge    = un.Age;
uAge2   = un.Age2;
uCons   = un.Cons;



%% Solve Worker's Problem

stagePay = (cMat + repmat(exp(w'),sizeV,1).*hMat)/rho+beta*PSI;
uCons = uCons + beta*PSI;

EMAX = zeros(sizeV,1);
EMAX_U = zeros(sizeV,1);

piPDF = zeros(sizeV,nOccs);
piPDF_U = zeros(sizeV,nOccs);


% 2. Terminal Age
D = exp(stagePay(sizeV - nTypes*nTen*nOccs + 1:end,:) - beta*PSI);
EMAX(sizeV - nTypes*nTen*nOccs + 1:end,1) = log(exp(uCons(state(sizeV - nTypes*nTen*nOccs + 1:end,2)) + uAge*nAge + uAge2*nAge^2 - beta*PSI) + sum(D,2));
piPDF(sizeV - nTypes*nTen*nOccs + 1:end,:) = bsxfun(@times,D,1./exp(EMAX(sizeV - nTypes*nTen*nOccs + 1:end,1)));

D = exp(uCost + stagePay(sizeV - nTypes*nTen*nOccs + 1:end,:) - beta*PSI);
EMAX_U(sizeV - nTypes*nTen*nOccs + 1:end,1) = log(exp(uCons(state(sizeV - nTypes*nTen*nOccs + 1:end,2)) + uAge*nAge + uAge2*nAge^2 - beta*PSI) + sum(D,2));
piPDF_U(sizeV - nTypes*nTen*nOccs + 1:end,:) = bsxfun(@times,D,1./exp(EMAX_U(sizeV - nTypes*nTen*nOccs + 1:end,1)));


% 3. Preceding Ages
for a = 2:nAge
    beg = sizeV - a*nTypes*nTen*nOccs + 1;
    fin = sizeV - (a-1)*nTypes*nTen*nOccs;
    
    E = EMAX(beg+nTypes*nTen*nOccs:fin+nTypes*nTen*nOccs);        % Continuation values for employment
    E_U = EMAX_U(beg+nTypes*nTen*nOccs:fin+nTypes*nTen*nOccs);    % Continuation values for unemployment
    D = exp(stagePay(beg:fin,:) + beta*E(T));                     % Numerators
    
    EMAX(beg:fin) = log(exp(uCons(state(beg:fin,2)) + uAge*(nAge-a+1) + uAge2*(nAge-a+1)^2 + beta*E_U) + sum(D,2));         % New Continuation Values
    piPDF(beg:fin,:) = bsxfun(@times,D,1./exp(EMAX(beg:fin,1)));             % Probability of move
    
    D = exp(uCost + stagePay(beg:fin,:) + beta*E(T));                        % Numerators if unemployed
    
    EMAX_U(beg:fin) = log(exp(uCons(state(beg:fin,2)) + uAge*(nAge-a+1) + uAge2*(nAge-a+1)^2 + beta*E_U(T_U)) + sum(D,2));  % New Unemp. Continuation Values
    piPDF_U(beg:fin,:) = bsxfun(@times,D,1./exp(EMAX_U(beg:fin,1)));
end

%% Go through histogram... [Can this be vectorized?]
% 2 histograms -- workers and then updating workers for next period
HS = zeros(nOccs,1);

for o = 1:nOccs
    A = dot(L0.*exp(bAge(o)*state(:,1) + bAge2(o)*(state(:,1).^2) + bTen(o)*(state(:,3).*(state(:,4)==o)) + bType(o,state(:,2))' + bRMSE(o)^2/2),piPDF(:,o));
    B = dot(LU0.*exp(bAge(o)*state(:,1) + bAge2(o)*(state(:,1).^2) + bTen(o)*(state(:,3).*(state(:,4)==o)) + bType(o,state(:,2))' + bRMSE(o)^2/2),piPDF_U(:,o));
    HS(o) = A+B;
end
