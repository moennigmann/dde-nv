classdef DDENLP < handle
    %this class is intended to for the robust stable steady state o
    %   Detailed explanation goes here
    
    % NV-status:
    % 1: input initial guess and might not solve the equations
    % 2: parameter is on critical manifold
    % 3: parameter is closest critical point
    % 4: Normal Vector at closest critical point was found
    % 5: connection of closest critical point and nominal point by normal
    % vector verified
    % 6: all nonlinear constraints have been concatenated
    % >6: values were handed back by optimization algorithm
    
    properties
        aCostFunction % cost function of the main steady state optimization problem
        problemDDE % object of class DDE, contains all relevant DDE information
        
        vars % structure containing the
        % cell array, find entries with  ind=find(ismember...
        
        stStCon = []; % of class equality constraints
        
        lowerBoxCons % low boundaries of box constraints
        upperBoxCons % upper boundaries of box contraints
        minDist = 1; % parameter uncertainy
        maxAllowedRealPart = 0; % the highest eigenvalue real part which is allowed.
        NVCon = []; % of class NVCon
        
        allNLEqConstraints % collects all nonlinear equality constraints
        allNLIneqConstraints = @(y)[]; % collects all nonlinear inequality constraints
        
        fixedUncertParamIndex % contains the indices of parameters that are uncertain, but fixed
        
        nlcon % is filled with the nonlinear constraints (inequality and equality)
        
        % numerical properties
        optionsInitEqCons = optimoptions('fsolve','Algorithm','levenberg-marquardt');
        optionsInitOptim = optimoptions('fmincon','Algorithm','active-set');
        optionsMainOptim = optimoptions('fmincon','Algorithm','active-set');
        
        % properties for stability/robustness check
        numMinEig=[];
        verifyStabPoints = 2; %  $$ {verifyStabPoints}^{nAlpha} $$ random points are generated and evaluated
        useLHS=1; % use latin hypercube sample for evaluating the uncertainty region
    end
    
    properties(SetAccess=protected)
        nX % number of states
        nAlpha % number of uncertain parameters
        nP = 0 % number of optimization parameters without uncertainty
        occupiedVars = 0; % number of variables in whole optimization problem
        occupiedEqs = 0; % number of equality constraints in whole optimization problem
        occupiedIneqs = 0; % number of inequality constraints  in whole optimization problem
        initVal % initial value, length(initVal)==occupiedVars
        optimVal % optimal value, should be a vector of  length(optimVal)==occupiedVars
        optJ % cost at optimum
        exitflag % exitflag of main optimization
        optimOutput % status message of main optimization
        lambda % lagrange multiplier of main optimization
        status = 0 % 0: nothing initialized
    end
    
    methods
        function aDDENLP = DDENLP(aCostFunction,aDDE,xNomGuess,stateLB,stateUB,uncertParamLB,uncertParamUB,varargin) % constructor
            % CONSTRUCTOR OF THE CLASS DDENLP
            % the simplest properties are assigned:
            
            
            % make sure all inputs are defined
            
            certOptParamLB=[];
            certOptParamUB=[];
            
            if nargin>=9
                certOptParamLB=varargin{1};
                certOptParamUB=varargin{2};
            end
            
            % the basics: cost function and dynamics:
            aDDENLP.aCostFunction = aCostFunction;
            aDDENLP.problemDDE = aDDE;
            
            %% box constraints
            % check if they are defined. If not: use default box
            % constraints
            
            % state
            if all(size(stateLB)==[xNomGuess.nVar,1])
                aDDENLP.lowerBoxCons = stateLB;
            else
                warning('lower boundaries for states have to be handed over as a column vector with as many entries as states in the DDE. Using a vector of -Inf instead')
                aDDENLP.lowerBoxCons = -Inf(xNomGuess.nVar,1);
            end
            
            if all(size(stateUB)==[xNomGuess.nVar,1])
                aDDENLP.upperBoxCons = stateUB;
            else
                warning('upper boundaries for states have to be handed over as a column vector with as many entries as states in the DDE. Using a vector of Inf instead')
                aDDENLP.upperBoxCons = Inf(xNomGuess.nVar,1);
            end
            
            % uncertain parameters
            if all(size(uncertParamLB)==[aDDE.uncParam.nVar,1])
                aDDENLP.lowerBoxCons = [aDDENLP.lowerBoxCons;uncertParamLB];
            else
                warning('lower boundaries for uncertain parameters have to be handed over as a column vector with as many entries as parameters in the DDE. Using a vector of -Inf instead')
                aDDENLP.lowerBoxCons = [aDDENLP.lowerBoxCons;-Inf(aDDE.uncParam.nVar,1)];
            end
            
            if all(size(uncertParamUB)==[aDDE.uncParam.nVar,1])
                aDDENLP.upperBoxCons = [aDDENLP.upperBoxCons;uncertParamUB];
            else
                warning('upper boundaries for uncertain parameters have to be handed over as a column vector with as many entries as parameters in the DDE. Using a vector of Inf instead')
                aDDENLP.upperBoxCons = [aDDENLP.upperBoxCons;Inf(aDDE.uncParam.nVar,1)];
            end
            
            % certain optimization parameters
            if all(size(certOptParamLB)==[aDDE.certOptParam.nVar,1])
                aDDENLP.lowerBoxCons = [aDDENLP.lowerBoxCons;certOptParamLB];
            else
                warning('lower boundaries for certain optimization parameters have to be handed over as a column vector with as many entries as parameters in the DDE. Using a vector of -Inf instead')
                aDDENLP.lowerBoxCons = [aDDENLP.lowerBoxCons;-Inf(aDDE.certOptParam.nVar,1)];
            end
            
            if all(size(certOptParamUB)==[aDDE.certOptParam.nVar,1])
                aDDENLP.upperBoxCons = [aDDENLP.upperBoxCons;certOptParamUB];
            else
                warning('upper boundaries for certain optimization parameters have to be handed over as a column vector with as many entries as parameters in the DDE. Using a vector of Inf instead')
                aDDENLP.upperBoxCons = [aDDENLP.upperBoxCons;Inf(aDDE.certOptParam.nVar,1)];
            end
            
            
            
            %% variables
            % the number of variables and states (the number of delays is
            % saved in the DDE object too avoid unnecessary redundance)
            aDDENLP.nX = xNomGuess.nVar;
            aDDENLP.nAlpha = aDDE.uncParam.nVar;
            aDDENLP.nP = aDDE.certOptParam.nVar;
            
            % guess for steady state
            aDDENLP.vars.nominal.x = xNomGuess;
            
            % nominal parameters
            aDDENLP.vars.nominal.alpha = aDDE.uncParam;
            
            aDDENLP.vars.nominal.p = aDDE.certOptParam;
            
            % we start without critical variables (they are added together
            % with NV-Constraints)
            aDDENLP.vars.critical = [];
            
            %% steady state constraint
            % add steady state constraint
            aDDENLP.stStCon = StStConstraint(aDDENLP.problemDDE,aDDENLP.vars.nominal);
            aDDENLP.initVal = [aDDENLP.stStCon.vars.x.values; aDDENLP.stStCon.vars.alpha.values];
            aDDENLP.occupiedVars = aDDENLP.occupiedVars+aDDENLP.stStCon.vars.x.nVar+aDDENLP.stStCon.vars.alpha.nVar+aDDENLP.vars.nominal.p.nVar;
            aDDENLP.occupiedEqs = aDDENLP.vars.nominal.x.nVar;
            aDDENLP.allNLEqConstraints = [aDDENLP.allNLEqConstraints; aDDENLP.stStCon.conFun]; %
            
        end
        
        function initializeStSt(aDDENLP)
            % the initial guess for the steady state is corrected to find
            % an actual steady state
            options=aDDENLP.optionsInitEqCons;
            aDDENLP.stStCon = initStStConstraint(aDDENLP.stStCon,options);
        end
        
        function addNVCon(aDDENLP,type,augSysHandle,nVSysHandle,xGuess,alphaGuess,pGuess)
            %point type (stst,fold,modfold,hopf,modhopf)
            % exist stst existAug existNV existConnection
            % entries
            
            nVvars = varCollection( type, aDDENLP.occupiedVars, xGuess, alphaGuess, pGuess);
            
            newNVCon = NVConstraint(aDDENLP,type,augSysHandle,nVSysHandle,nVvars);
            aDDENLP.NVCon = [aDDENLP.NVCon; newNVCon];
            
            aDDENLP.occupiedVars = aDDENLP.occupiedVars+...
                +newNVCon.nVarAugSys...
                +newNVCon.nVarNVSys...
                +1;
            
            aDDENLP.occupiedEqs = aDDENLP.occupiedEqs+...
                +newNVCon.nEqs;
            %             newNVCon.nEqs
            %             size(augSysHandle(zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1)))
            %             size(nVSysHandle(zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1),zeros(200,1)))
            %
            aDDENLP.occupiedIneqs = aDDENLP.occupiedIneqs+...
                + 1;
            
            aDDENLP.vars.critical = [aDDENLP.vars.critical;newNVCon.vars];
            
            aDDENLP.status=min(aDDENLP.status,newNVCon.status);
            
            
        end
        
        function initNVCons(aDDENLP)
            % this function conducts the steps necessary for the
            % initialization of all normal vector constraints i.e.
            % - finding a critical point
            % - finding the closesest critical point
            % - check if the closest critical point exists allready in
            % another normal vector constraint
            % - finding the normal vector at the closest critical point
            % - connect closest critical point and nominal point
            
            alphaNom=aDDENLP.vars.nominal.alpha;
            
            for i=1:length(aDDENLP.NVCon)
                aDDENLP.NVCon(i).findManifoldPoint(aDDENLP.NVCon(i).vars);
                aDDENLP.NVCon(i).findClosestCriticalPoint(alphaNom);
                for j=1:i-1
                    % if another closest critical point is very close (i.e.
                    % possibly identical, a warning is shown
                    if norm(aDDENLP.NVCon(i).vars.alpha.values-aDDENLP.NVCon(j).vars.alpha.values,2)<1e-3
                        warning('closest critical points of NVCon(%d) and NVCon(%d) conincide!\n',j,i)
                        return
                    end
                end
                aDDENLP.NVCon(i).findNormalVector(alphaNom);
                aDDENLP.NVCon(i).findConnection(alphaNom);
            end
            
            
        end
        
        function evaluateStatus(aDDENLP)
            % checks the status of all existing constraints and assign the
            % status of the whole optimization problem
            if aDDENLP.stStCon.status==1
                aDDENLP.status=min(arrayfun(@(x)x.status,aDDENLP.NVCon));
            end
        end
        
        function concatConstraints(aDDENLP)
            % concatenates the existing constraints
            % fmincon expects nonlinear constraints to have the form
            % [c,ceq]=nlcon(x). The different nonlinear constraints are
            % still distributed over various objects. They are all
            % collected and a combined constraint function is formed
            
            %start with steaty state constraint
            nVEq=@(y)aDDENLP.stStCon.conFun(y)';
            nVIneq=@(y)[];
            
            % and add attach the normal vector constraints
            for i=1:size(aDDENLP.NVCon,1)
                
                nVEq =@(y) [...
                    nVEq(y);
                    aDDENLP.NVCon(i).conFun(y)];
                
                nVIneq=@(y)[...
                    nVIneq(y);
                    aDDENLP.NVCon(i).inequalities(y)];
                
            end
            
            % save the constraints as object properties
            aDDENLP.allNLEqConstraints=nVEq;
            aDDENLP.allNLIneqConstraints=nVIneq;
            
            % put equality and inequality constraints together
            aDDENLP.nlcon=@(x)deal(aDDENLP.allNLIneqConstraints(x),aDDENLP.allNLEqConstraints(x));
            
        end
        
        function concatInitPoints(aDDENLP)
            % the critical points and normal vectors have been stored
            % distributed in different object. This function gathers all of
            % them and saves them as a variable to be used in the main
            % optimization
            
            aDDENLP.initVal=[...
                aDDENLP.stStCon.vars.x.values;...
                aDDENLP.stStCon.vars.alpha.values;...
                aDDENLP.vars.nominal.p.values];
            
            for i=1:length(aDDENLP.NVCon)
                
                newVars = [...
                    aDDENLP.NVCon(i).vars.x.values;...
                    aDDENLP.NVCon(i).vars.alpha.values;...
                    aDDENLP.NVCon(i).vars.omega.values;...
                    aDDENLP.NVCon(i).vars.w1.values;...
                    aDDENLP.NVCon(i).vars.w2.values;...
                    aDDENLP.NVCon(i).vars.v1.values;...
                    aDDENLP.NVCon(i).vars.v2.values;...
                    aDDENLP.NVCon(i).vars.g1.values;...
                    aDDENLP.NVCon(i).vars.g2.values;...
                    aDDENLP.NVCon(i).vars.u.values;...
                    aDDENLP.NVCon(i).vars.r.values;...
                    aDDENLP.NVCon(i).vars.l.values];
                aDDENLP.initVal=[aDDENLP.initVal;newVars];
                
            end
        end
        
        function checkConstraints(aDDENLP)
            % displays which constraints might be violated at initial point
            %
            %             fprintf('\nLower Box Contstraints: \n')
            %             size(aDDENLP.lowerBoxCons)
            %             size( aDDENLP.initVal)
            %             disp(aDDENLP.lowerBoxCons < aDDENLP.initVal) % low boundaries of box constraints
            %             fprintf('\nLower Box Contstraints: \n')
            %             disp(aDDENLP.upperBoxCons < aDDENLP.initVal) % upper boundaries of box contraints
            %
            fprintf('\n nonlinear inequality constraints: \n')
            disp(aDDENLP.allNLIneqConstraints(aDDENLP.initVal)<0)
        end
        
        function runOptim(aDDENLP,varargin)
            % this function starts the main optimization without
            % intermediate stability checks
            
            if nargin > 1
                options = varargin{1};
            else
                options = aDDENLP.optionsMainOptim;
            end
            
            % store, how this function was called (for later use)
            callerFunction = dbstack;
            
            % check status
            aDDENLP.evaluateStatus();
            
            if aDDENLP.status<5
                error('constraints appear to be not completely initialized')
            end
            
            % define shorter names for the optimization arguments
            J = aDDENLP.aCostFunction;
            x0 = aDDENLP.initVal;
            lb = aDDENLP.lowerBoxCons;
            ub = aDDENLP.upperBoxCons;
            
            nonlinCon = aDDENLP.nlcon;
            
            % prepare constraints
            lb = [lb; -Inf(length(x0)-length(lb),1)];
            
            ub = [ub; Inf(length(x0)-length(ub),1)];
            
            if isempty(aDDENLP.nlcon)
                warning('no nonlinear constraints found. Did you concatenate the nonlinear constraints with the function ''aDDENLP = concatConstraints(aDDENLP)''?')
            end
            
            Aeq=eye(length(x0));
            Aeq=Aeq(aDDENLP.fixedUncertParamIndex,:);
            beq=x0(aDDENLP.fixedUncertParamIndex);
            
            [c,ceq]=nonlinCon(x0);
            
            if any(c>0)
                if ~strcmp(callerFunction(2).name, 'DDENLP.runOptimWithStabChecks')
                    warning('nonlinear inequality constraints violated at initial point')
                end
            end
            
            %             if (norm(ceq,2)>options.TolCon)
            %                 fprintf('%d\n',ceq);
            %                 error('nonlinear equality constraints violated at initial point (measured by TolCon = %d).',(options.TolCon))
            %             end
            
            % run optimization
            [aDDENLP.optimVal,aDDENLP.optJ,aDDENLP.exitflag,aDDENLP.optimOutput,aDDENLP.lambda] = ...
                fmincon(J,x0,[],[],Aeq,beq,lb,ub,nonlinCon,options);
            
            aDDENLP.vars.nominal.x.values = aDDENLP.optimVal(aDDENLP.vars.nominal.x.index);
            aDDENLP.vars.nominal.alpha.values = aDDENLP.optimVal(aDDENLP.vars.nominal.alpha.index);
            aDDENLP.vars.nominal.p.values = aDDENLP.optimVal(aDDENLP.vars.nominal.p.index);
            
            if aDDENLP.exitflag>0
                aDDENLP.status = 6;
                fprintf('\nmain optimization finished!\n');
            else
                if ~strcmp(callerFunction(2).name, 'DDENLP.runOptimWithStabChecks')
                    fprintf('\nmain optimization not successful, exitflag was %d \n', aDDENLP.exitflag);
                end
            end
        end
        
        
        %% WORK IN PROGRESS: STARTING MARK
        
        function [init, final, maxEig] = runOptimWithStabChecks(aDDENLP,nIterBetweenStabChecks)
            % this function starts the main optimization with
            % intermediate stability checks
            
            options = aDDENLP.optionsMainOptim;
            maxIter = options.MaxIter;
            
            options.MaxIter = nIterBetweenStabChecks;
            options.Display = 'off';
            
            for ii = 1:ceil(maxIter/nIterBetweenStabChecks)
                aDDENLP.concatInitPoints();
                aDDENLP.runOptim(options);
                init = aDDENLP.deconstructInit()
                final = aDDENLP.deconstructOptimum();
                
                [maxEig,~] = aDDENLP.checkStabilityPoint('nominal');
                
                if real(maxEig) > aDDENLP.maxAllowedRealPart
                    warning('System lost desired stability properties during optimization. You can use the properies DDENLP.initVal and DDENLP.opitimVal to find out were stability was lost')
                    break
                end
                % leave loop if optimization finished
                if aDDENLP.exitflag > 0
                    break
                end
                
            end
        end
        
        function runOptimAddingNewManifolds(aDDENLP,nIterBetweenStabChecks)
            
            
            for ii = 1:10
                [init,final,maxEig] = aDDENLP.runOptimWithStabChecks(nIterBetweenStabChecks);
                init.nominal.alpha.values
                final.nominal.alpha.values
                if real(maxEig) > aDDENLP.maxAllowedRealPart
                    if imag(maxEig) == 0
                        subtype = 'fold';
                    else
                        subtype = 'hopf';
                    end
                    if aDDENLP.maxAllowedRealPart ~= 0
                        type = ['mod',subtype];
                    else
                        type = subtype;
                    end
                    disp(type)
                    [handleManifold,handleNV] = aDDENLP.problemDDE.getHandles(type);
                    intermediatePoint = aDDENLP.findManifoldPointOnLine(type, handleManifold, init, final);
                    
                    %add new NVCON                    
                    aDDENLP.addNV(type,handleManifold,handleNV,intermediatePoint.x,intermediatePoint.alpha,init.p)
                    % and initialize it:
                    aDDENLP.NVCon(end).findManifoldPoint(aDDENLP.NVCon(end).vars);
                    aDDENLP.NVCon(end).findClosestCriticalPoint(alphaNom);
                    aDDENLP.NVCon(end).findNormalVector(alphaNom);
                    aDDENLP.NVCon(end).findConnection(alphaNom);
                else
                    break
                end
            end
        end
        
        
        function  intermediatePoint = findManifoldPointOnLine(aDDENLP, type, manifoldHandle, point1, point2 )
            % This method looks for a crossing of a critical manifold between two given
            % points
            
            alphaInterp = @(lambda)point1.alpha.values-lambda*(point1.alpha.values-point2.alpha.values);
            pInterp = @(lambda)point1.p.values-lambda*(point1.p.values-point2.p.values);
            
            switch type
                case 'fold'
                    funHandle = @(x,lambda,w1)manifoldHandle(x,alphaInterp(lambda),pInterp(lambda),w1);
                    rhs = @(y)funHandle(y(1:point1.x.nVar),...
                        y(point1.x.nVar+1),...
                        y(point1.x.nVar+2:2*point1.x.nVar+1));
                case 'modfold'
                    funHandle = @(x,lambda,w1)manifoldHandle(x,alphaInterp(lambda),pInterp(lambda),w1);
                    rhs = @(y)funHandle(y(1:point1.x.nVar),...
                        y(point1.x.nVar+1),...
                        y(point1.x.nVar+2:2*point1.x.nVar+1));
                case 'hopf'
                    funHandle = @(x,lambda,omega,w1,w2)manifoldHandle(x,alphaInterp(lambda),pInterp(lambda),omega,w1,w2);
                    rhs = @(y)funHandle(y(1:point1.x.nVar),...
                        y(point1.x.nVar+1),...
                        y(point1.x.nVar+2),...
                        y(point1.x.nVar+3:2*point1.x.nVar+2),...
                        y(2*point1.x.nVar+3:3*point1.x.nVar+2));
                case 'modhopf'
                    funHandle = @(x,lambda,omega,w1,w2)manifoldHandle(x,alphaInterp(lambda),pInterp(lambda),omega,w1,w2);
                    rhs = @(y)funHandle(y(1:point1.x.nVar),...
                        y(point1.x.nVar+1),...
                        y(point1.x.nVar+2),...
                        y(point1.x.nVar+3:2*point1.x.nVar+2),...
                        y(2*point1.x.nVar+3:3*point1.x.nVar+2));
                otherwise
                    error('unknown type requested')
            end
            
            
            y0 = [ point1.x.values-0.5*(point1.x.values-point2.x.values);...
                0.5;...
                0;...
                ones(size(point1.x.values));...
                zeros(size(point1.x.values)) ];
            
            
            [y,~,exitflag] = fsolve(rhs,y0,aDDENLP.optionsInitEqCons);
            
            if exitflag > 0
                x = y(1:point1.x.nVar);
                lambda = y(point1.x.nVar+1);
                alpha = alphaInterp( lambda );
                p =  pInterp( lambda);
                omega = y(point1.x.nVar+2);
                w1 = y(point1.x.nVar+3:2*point1.x.nVar+2);
                w2 = y(2*point1.x.nVar+3:3*point1.x.nVar+2);
                
                intermediatePoint.x = VariableVector(x,Inf,point1.x.names);
                intermediatePoint.alpha = VariableVector(alpha, Inf,point1.alpha.names);
                intermediatePoint.p = VariableVector(p, Inf,point1.p.names);
                intermediatePoint.omega = VariableVector(omega, Inf);
                intermediatePoint.w1 = VariableVector(w1, Inf);
                intermediatePoint.w2 = VariableVector(w2, Inf);
            else
                warning('did not find critical point on the given line, fsolve exitflag was %d',exitflag)
            end
            
        end
        
        
        %% FROM HERE UPWARDS: ENDING MARK
        
        function [varargout] = deconstructInit(aDDENLP)
            if nargout > 0
                varargout{1} = 1;
            end
            init = aDDENLP.initVal;
            outputVars = aDDENLP.vars;
            outputVars.nominal.x.values = init(outputVars.nominal.x.index);
            outputVars.nominal.alpha.values = init(outputVars.nominal.alpha.index);
            outputVars.nominal.p.values = init(outputVars.nominal.p.index);
            
            for i=1:length(aDDENLP.NVCon)
                outputVars.critical(i).x.values = init(outputVars.critical(i).x.index);
                outputVars.critical(i).alpha.values = init(outputVars.critical(i).alpha.index);
                outputVars.critical(i).omega.values = init(outputVars.critical(i).omega.index);
                outputVars.critical(i).w1.values = init(outputVars.critical(i).w1.index);
                outputVars.critical(i).w2.values = init(outputVars.critical(i).w2.index);
                outputVars.critical(i).v1.values = init(outputVars.critical(i).v1.index);
                outputVars.critical(i).v2.values = init(outputVars.critical(i).v2.index);
                outputVars.critical(i).g1.values = init(outputVars.critical(i).g1.index);
                outputVars.critical(i).g2.values = init(outputVars.critical(i).g2.index);
                outputVars.critical(i).u.values = init(outputVars.critical(i).u.index);
                outputVars.critical(i).r.values = init(outputVars.critical(i).r.index);
                outputVars.critical(i).l.values = init(outputVars.critical(i).l.index);
            end
            aDDENLP.vars=outputVars;
        end
        
        function varargout = deconstructOptimum(aDDENLP)
            
            if nargout > 0
                varargout{1} = 1;
            end

            opt = aDDENLP.optimVal;
            outputVars = aDDENLP.vars;
            outputVars.nominal.x.values = opt(outputVars.nominal.x.index);
            outputVars.nominal.alpha.values = opt(outputVars.nominal.alpha.index);
            outputVars.nominal.p.values = opt(outputVars.nominal.p.index);
            
            for i=1:length(aDDENLP.NVCon)
                outputVars.critical(i).x.values = opt(outputVars.critical(i).x.index);
                outputVars.critical(i).alpha.values = opt(outputVars.critical(i).alpha.index);
                outputVars.critical(i).omega.values = opt(outputVars.critical(i).omega.index);
                outputVars.critical(i).w1.values = opt(outputVars.critical(i).w1.index);
                outputVars.critical(i).w2.values = opt(outputVars.critical(i).w2.index);
                outputVars.critical(i).v1.values = opt(outputVars.critical(i).v1.index);
                outputVars.critical(i).v2.values = opt(outputVars.critical(i).v2.index);
                outputVars.critical(i).g1.values = opt(outputVars.critical(i).g1.index);
                outputVars.critical(i).g2.values = opt(outputVars.critical(i).g2.index);
                outputVars.critical(i).u.values = opt(outputVars.critical(i).u.index);
                outputVars.critical(i).r.values = opt(outputVars.critical(i).r.index);
                outputVars.critical(i).l.values = opt(outputVars.critical(i).l.index);
            end
            aDDENLP.vars=outputVars;
            
        end
        
        function [maxRealPart,eigs] = checkStabilityPoint(aDDENLP,point)
            if ischar(point)
                type=point;
                switch point
                    case 'nominal'
                        point=aDDENLP.vars.nominal;
                    case 'critical'
                        point=aDDENLP.vars.critical;
                    case 'vertex'
                        point=aDDENLP.vars.vertex;
                    case 'random'
                        point=aDDENLP.vars.random;
                    otherwise
                        error('unknown point type requested')
                end
            end
            maxRealPart = NaN(size(point,1),1);
            
            p=aDDENLP.vars.nominal.p.values;
            
            myDDE = @(xx,alpha)aDDENLP.problemDDE.rhs(xx(:,1),xx(:,2:end),alpha,p)';
            myNtau = aDDENLP.problemDDE.ntau;
            myDelays = @(k,xx,alpha)(k==(1:myNtau))*aDDENLP.problemDDE.delays(xx,alpha,p);
            
            funcs = set_funcs('sys_rhs',myDDE,...
                'sys_ntau',myNtau,...
                'sys_tau',myDelays);
            
            for ii = 1 : size(point,1);
                [maxRealPartTMP, xCorr, eigs] = checkStability(funcs,point(ii).alpha.values,point(ii).x.values,aDDENLP.numMinEig);
                if isempty(maxRealPartTMP)
                    error('did not find an eigevalue with Re(lambda)>aDDENLP.numMinEig. Please chose a smaller value for aDDENLP.numMinEig')
                end
                maxRealPart(ii) = maxRealPartTMP;
                switch type
                    case 'nominal'
                        aDDENLP.vars.nominal(ii,1).maxEig = maxRealPart(ii);
                    case 'critical'
                        aDDENLP.vars.critical(ii,1).maxEig = maxRealPart(ii);
                    case 'vertex'
                        aDDENLP.vars.vertex(ii,1).x.values = xCorr;
                        aDDENLP.vars.vertex(ii,1).maxEig = maxRealPart(ii);
                    case 'random'
                        aDDENLP.vars.random(ii,1).x.values = xCorr;
                        aDDENLP.vars.random(ii,1).maxEig = maxRealPart(ii);
                end
            end
            %             switch type
            %                 case 'nominal'
            %                     aDDENLP.vars.nominal(ii,1).maxEig=maxRealPart(ii);
            %                 case 'critical'
            %                     aDDENLP.vars.critical(ii,1).maxEig=maxRealPart(ii);
            %                 case 'vertex'
            %                     aDDENLP.vars.vertex(ii,1).x.values=xCorr;
            %                     aDDENLP.vars.vertex(ii,1).maxEig=maxRealPart(ii);
            %                 case 'random'
            %                     aDDENLP.vars.random(ii,1).x.values=xCorr;
            %                     aDDENLP.vars.random(ii,1).maxEig=maxRealPart(ii);
            %             end
        end
        
        function maxRealPart = checkStabilityAtVertices(aDDENLP)
            
            myNAlpha=aDDENLP.nAlpha;
            
            vertices(2^myNAlpha,1)=struct();
            
            errors = dec2bin(0:2^myNAlpha-1);
            errorvector = -ones(2^myNAlpha,myNAlpha);
            
            for i=1:myNAlpha
                errorvector(:,i)=errorvector(:,i)+2*bin2dec(errors(:,i));
            end
            
            
            for i=1:2^myNAlpha
                vertices(i).alpha.values = aDDENLP.vars.nominal.alpha.values+errorvector(i,:)';
                vertices(i).x.values = aDDENLP.vars.nominal.x.values;
                
                aDDENLP.vars.vertex(i,1).alpha=VariableVector(vertices(i).alpha.values,NaN,{'alphaAtVertex'});
                aDDENLP.vars.vertex(i,1).x=VariableVector(vertices(i).x.values,NaN,{'xAtVertex'});
            end
            
            maxRealPart = checkStabilityPoint(aDDENLP,'vertex');
        end
        
        function maxRealPart = checkStabilityAtRandom(aDDENLP, baseForNumberOfPoints, seedForRandomNumbers, varargin)
            % calculates max(real(eigenvalues)) for random points within
            % the uncertainty region
            
            if nargin > 2
                rng(varargin{2});
            else
                if nargin == 1
                    baseForNumberOfPoints = 2;
                end
            end
            
            nPoints = ceil(baseForNumberOfPoints^aDDENLP.nAlpha);
            
            randPoint(nPoints,1)=struct();
            
            if aDDENLP.useLHS
                pointVector = 2*lhsdesign(nPoints,aDDENLP.nAlpha)-ones(nPoints,aDDENLP.nAlpha);
            else
                pointVector = 2*rand(nPoints,aDDENLP.nAlpha)-ones(nPoints,aDDENLP.nAlpha);
            end
            
            for i=1:nPoints
                randPoint(i).alpha.values = aDDENLP.vars.nominal.alpha.values+pointVector(i,:)';
                randPoint(i).x.values = aDDENLP.vars.nominal.x.values;
                
                aDDENLP.vars.random(i,1).alpha=VariableVector(randPoint(i).alpha.values,NaN,{'alphaAtVertex'});
                aDDENLP.vars.random(i,1).x=VariableVector(randPoint(i).x.values,NaN,{'xAtVertex'});
            end
            
            [maxRealPart,~] = checkStabilityPoint(aDDENLP,'random');
        end
        
        function sol = ddesd(aDDENLP,point,history,tspan,options)
            
            myDDE=@(~,x,xtau)aDDENLP.problemDDE.rhs(x,xtau,point.alpha.values,point.p.values)';
            myDelays=@(t,x)t-aDDENLP.problemDDE.delays(x,point.alpha.values,point.p.values)';
            
            sol = ddesd(myDDE,myDelays,history,tspan,options) ;
        end
        
        
        %         generateDDE % f(t,x,xtau,alpha)
        %         generateDelayVec % tau(x,alpha)
        
        %         generateHopfMani
        %         generateHopfNV
        %
        %         generateFoldMani
        %         generateFoldNV
        %
        %         generateModHopfMani
        %         generateModHopfNV
        %
        %         generateModFoldMani
        %         generateModFoldNV
        
        %         verifyRobustness
        
        %         export2DDEBIFTOOL % nicht bei class DDE, damit Box
        %         constraints auch übergeben werden können
        
    end
    
end

