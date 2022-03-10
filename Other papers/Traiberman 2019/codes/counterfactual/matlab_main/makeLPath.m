function [LPath, LPathU, occDistPath, piCells, piUCells] = makeLPath(L0,LU0,realWage,VPath,VPathU,state,n,b,un,params,TAMAT,TAMAT_U,cMat,hMat,numT)


% Model Parameters
nTypes   =  n.Types;
nTen     =  n.Ten;
nAge     =  n.Age;
nOccs    =  n.Occs;
sizeV    =  nTypes*nTen*nOccs*nAge;

% Mincer Parameters
bAge  = b.age;
bAge2 = b.age2;
bTen  = b.Ten;
bType = b.Type;
bRMSE = b.RMSE;

LPath = zeros(sizeV,numT);
LPath(:,1) = L0;

LPathU = zeros(sizeV,numT);
LPathU(:,1) = LU0;

occDistPath = zeros(nOccs,numT);
occDistPath(:,1) = accumarray(state(:,4),L0.*exp(bRMSE(state(:,4)).^2/2 + bAge(state(:,4)).*state(:,1) + bAge2(state(:,4)).*(state(:,1).^2) + bTen(state(:,4)).*state(:,3) + bType(sub2ind(size(bType),state(:,4),state(:,2)))));

piCells = cell(numT-1,1);
piUCells = cell(numT-1,1);

for j = 2:numT;
        % Step 4: Calculate new distribution of workers
        [occDistPath(:,j),piPDF, piPDF_U] = laborSupplyF_pf_new(realWage(:,j),VPath(:,j-1),VPathU(:,j-1),L0,LU0,n,un,b,params,TAMAT,TAMAT_U,cMat,hMat,state);
        
        piCells{j-1} = piPDF;
        piUCells{j-1} = piPDF_U;
        
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
        
        LPath(:,j) = L1;
        LPathU(:,j) = LU1;
        
        L0 = L1;
        LU0 = LU1;
end