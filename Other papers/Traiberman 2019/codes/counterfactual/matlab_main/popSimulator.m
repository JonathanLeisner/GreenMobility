function [history, switches, total] = popSimulator(LCDF, piPath,piUPath, state, simN, seed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to simulate histories for a random sample of the population
% after the model has been solved.
% 
% Inputs [Required]:
% 1. LCDF       = initial CDF across agents
% 2. piPath     = PDF across states (employed)
% 3. piUPath    = PDF across states (unemployed)
% 4. state      = matrix of states
% 
% Inputs [Optional]: 
% 1. simN       = number of simulations [default = 10000] 
% 2. seed       = seed for number stream [default = 479]
% 
% 
% Outputs:
% 1. history    = cell array of agent histories
% 2. switches   = time series of total occupation switches
% 3. total      = total number of folks still breathing
% 
% 
% Written by Sharon Traiberman
% First version: 11/29/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% # of simulations excluded
if nargin == 4
    simN = 100000;
    seed = 479;
end

if isempty(simN)
	simN = 100000;
end

% seed excluded
if nargin == 5

    seed = 479;
end


% Generate Stream
s = RandStream('mt19937ar','Seed',seed);
simDraws = rand(s,simN,36);

% Some preallocation
history = cell(simN,1);
switches = zeros(35,1);
total = zeros(35,1);


tic
for simn = 1:simN

    e = simDraws(simn,1);

    initInd = find(LCDF>e,1);
    if isempty(initInd)
        initInd = 1;
    end

    % State = (age, type, tenure, occ)
    history{simn} = zeros(35,9);
    history{simn}(1,:) = [simn,1,initInd,state(initInd,:),1,0];

    curAge = state(initInd,1);
    curType = state(initInd,2);
    curTen = state(initInd,3);
    curOcc = state(initInd,4);
    curInd = initInd;

    curEmp = 1;

    for t = 1:35
        e = simDraws(simn,t+1);
        curSwitch = 0;
        % Stop at retirement
        if curAge>=35
            break
        end
        total(t) = total(t)+1;

        % If employed
        if curEmp==1        
            curCDF = cumsum(piPath{t}(curInd,:));
        else
            curCDF = cumsum(piUPath{t}(curInd,:));
        end
        newOcc = find(curCDF>e,1);

        % 1. Updating Age [Easy]
        curAge = curAge+1;

        % 2. Updating Employment Status, Tenure and Occupation
        % Check if entering non-employment
        if isempty(newOcc)

            % Reset tenure if too long unemployed
            if curEmp==0
                curTen = 1;
            end

            % Set employment to 0
            curEmp=0;
        else
            
            % Set employment to 1
            curEmp=1;

            % Check if a switch to update tenure
            if curOcc~=newOcc
                curTen = 1;
                switches(t) = switches(t)+1;
                curSwitch = 1;
            else
                curTen = min(curTen+1,6);
            end

            % New occupation
            curOcc = newOcc;
        end

        % New Position
        curInd = find((state(:,1)==curAge) & (state(:,2)==curType) & (state(:,3)==curTen) & (state(:,4)==curOcc));
        history{simn}(t+1,:) = [simn,t+1,curInd,state(curInd,:),curEmp,curSwitch];  
    end 
end 