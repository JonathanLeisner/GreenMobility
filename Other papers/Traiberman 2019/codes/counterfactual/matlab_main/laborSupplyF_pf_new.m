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

function [HS, piPDF, piPDF_U, EMAX] = laborSupplyF_pf_new(w,V,VU,L0,LU0,n,un,b,params,TA,TA_U,cMat,hMat,state)

%% Prelminary Constants
% Mathematical Constants
PSI = -psi(1);

% Model Parameters
nTypes   =  n.Types;
nTen     =  n.Ten;
nAge     =  n.Age;
nOccs    =  n.Occs;
sizeV    =  nTypes*nTen*nOccs*nAge;

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

%% Calculate Stage Payoffs
stagePay = (cMat + repmat(exp(w'),sizeV,1).*hMat)/rho+beta*PSI;
uCons = uCons + beta*PSI;

%% Moving continuation values up one on the life cycle...
Vtemp = zeros(sizeV,1);
Vtemp(1:nTypes*nTen*nOccs*(nAge-1)) = V(nTypes*nTen*nOccs+1:end);
Vtemp(sizeV - nTypes*nTen*nOccs + 1:end,1) = -PSI;

VtempU = zeros(sizeV,1);
VtempU(1:nTypes*nTen*nOccs*(nAge-1)) = VU(nTypes*nTen*nOccs+1:end);
VtempU(sizeV - nTypes*nTen*nOccs + 1:end,1) = -PSI;

%% Solving for transition matrix
D = exp(stagePay + beta*Vtemp(TA));
EMAX = exp(uCons(state(:,2)) + uAge*state(:,1) + uAge2*state(:,1).^2 + beta*VtempU) + sum(D,2);
piPDF = bsxfun(@times,D,1./EMAX);


DU = exp(uCost + stagePay + beta*Vtemp(TA));
EMAX_U = exp(uCons(state(:,2)) + uAge*state(:,1) + uAge2*state(:,1).^2 + beta*VtempU(TA_U)) + sum(DU,2);
piPDF_U = bsxfun(@times,DU,1./EMAX_U);


%% Go through histogram...
HS = zeros(nOccs,1);

for o = 1:nOccs
    A = dot(L0.*exp(bAge(o)*state(:,1) + bAge2(o)*(state(:,1).^2) + bTen(o)*(state(:,3).*(state(:,4)==o)) + bType(o,state(:,2))' + bRMSE(o)^2/2),piPDF(:,o));
    B = dot(LU0.*exp(bAge(o)*state(:,1) + bAge2(o)*state(:,1).^2 + bTen(o)*(state(:,3).*(state(:,4)==o)) + bType(o,state(:,2))' + bRMSE(o)^2/2),piPDF_U(:,o));
    HS(o) = A+B;
end
