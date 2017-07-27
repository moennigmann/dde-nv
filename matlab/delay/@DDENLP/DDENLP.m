%> @file DDENLP.m
%> @brief this class desribes the whole optimization problem for delayed systems
%> 
%> @mainpage Documentation of the Toolbox DDENLP for robust steady state optimization of delayed systems.
%>
%> @section use Using this toolbox
%> This toolbox aids with the optimization of delayed systems using normal vector constraints for robust stability properties.
%>
%> The code documentation can be found in this file.
%>
%> Additionally, there are some example files.
%>
%> @section funding Funding 
%> The development of this toolbox was funded by  Deutsche Forschungsgemeinschaft (grant MO 1086/13).
%>
%> @author Jonas Otten and Martin Moennigmann
%> @date 18 Jul 2017
%>
%>======================================================================
%> @brief this class desribes the whole optimization problem for delayed systems
%
%>  this class is intended to for the robust stable steady state of
%>  delayed systems
%> 
%>     status:
%>     1: input initial guess and might not solve the equations
%>     2: parameter is on critical manifold
%>     3: parameter is closest critical point
%>     4: Normal Vector at closest critical point was found
%>     5: connection of closest critical point and nominal point by normal
%>     vector verified
%>     6: all nonlinear constraints have been concatenated
%>     >6: values were handed back by optimization algorithm
%>
%> @author Jonas Otten
%> @date 18 Jul 2017 ======================================================================




classdef DDENLP < handle
    %this class is intended to for the robust stable steady state of
    %delayed systems
      
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
        %> cost function of the main steady state optimization problem
        aCostFunction 
        %> object of class DDE, contains all relevant DDE information
        problemDDE
        %> structure containing the cell array of variables, find entries with  ind=find(ismember...
        vars 
        
        %> steady state constraint class StStConstraint
        stStCon = []; 
        
        %> lower boundaries of box constraints
        lowerBoxCons
        %> lower boundaries of box constraints
        upperBoxCons
        %> parameter uncertainty
        minDist = 1;
        %> the highest eigenvalue real part which is allowed 
        maxAllowedRealPart = 0;
        %> accepted eigenvalues in the closed right halfplane (for eigenvalues on the imaginary axis)
        allowedEigsInClosedRightHP = 0;
        
        %> normal vector constraints, vector of objects of class NVConstraint
        NVCon = [];
        
        %> function handle that collects all nonlinear equality constraints
        allNLEqConstraints 
        %> function handle that collects all nonlinear inequality constraints
        allNLIneqConstraints = @(y)[]; 
        
        %> contains the indices of parameters that are uncertain, but fixed
        fixedUncertParamIndex 
        
        %> function handle with the nonlinear constraints (inequality and equality)
        nlcon 
        
        % numerical properties
        %> options for initialization numerics (fsolve)
        optionsInitEqCons = optimoptions('fsolve','Algorithm','levenberg-marquardt');
        %> options for initialization numerics (auxiliary optimization)
        optionsInitOptim = optimoptions('fmincon','Algorithm','active-set');
        %> options for the main optimization
        optionsMainOptim = optimoptions('fmincon','Algorithm','active-set');
        
        % properties for stability/robustness check
        %> real number, for DDE-BIFTOOL numerics
        numMinEig=[];
        %>  $$ {verifyStabPoints}^{nAlpha} $$ random points are generated and evaluated
        verifyStabPoints = 2; 
        %> flag for using latin hypercube sample for evaluating the uncertainty region
        useLHS=1;
    end

    properties(SetAccess=protected)
        %> number of states
        nX
        %> number of uncertain parameters 
        nAlpha
        %> number of optimization parameters without uncertainty
        nP = 0 
        %> number of variables in whole optimization problem 
        occupiedVars = 0;
        %> number of equality constraints in whole optimization problem
        occupiedEqs = 0; 
        %> number of inequality constraints  in whole optimization problem
        occupiedIneqs = 0; 
        %> initial value, length(initVal)==occupiedVars
        initVal
        %> optimal value, should be a vector of  length(optimVal)==occupiedVars
        optimVal
        %> cost at optimum
        optJ
        %> exitflag of main optimization
        exitflag
        %> status message of main optimization
        optimOutput
        %> lagrange multiplier of main optimization
        lambda
        %> status of this optimization problem, (key in main description)
        status = 0 
    end
    
    methods
        
    % ======================================================================
    %> @brief Class constructor
    %>
    %> This function constructs instances of the class DDENLP
    %>
    %> @param aCostFunction function handle for cost function
    %> @param aDDE underlying DDE, instance of class DDE
    %> @param xNomGuess collection of instances of VariableVector, nominal
    %>        state guess
    %> @param stateLB lower bounds for states
    %> @param stateUB upper bounds for states
    %> @param uncertParamLB lower bounds for uncertain optimization parameters
    %> @param uncertParamUB upper bounds for uncertain optimization parameters
    %> @param certOptParamLB lower bounds for certain optimization parameters (optional input)
    %> @param certOptParamUB upper bounds for certain optimization parameters (optional input)
    % 
    %> @return instance of the DDENLP class.
    % ======================================================================
        
        function aDDENLP = DDENLP( aCostFunction, aDDE, xNomGuess, stateLB, stateUB, uncertParamLB, uncertParamUB, certOptParamLB, certOptParamUB, varargin ) % constructor
            % the simplest properties are assigned:
            
            
            % make sure all inputs are defined
            if nargin<9
                certOptParamLB=[];
                certOptParamUB=[];
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
                if aDDE.certOptParam.nVar > 0
                    warning('lower boundaries for certain optimization parameters have to be handed over as a column vector with as many entries as parameters in the DDE. Using a vector of -Inf instead')
                end
                aDDENLP.lowerBoxCons = [aDDENLP.lowerBoxCons;-Inf(aDDE.certOptParam.nVar,1)];
            end
            
            if all(size(certOptParamUB)==[aDDE.certOptParam.nVar,1])
                aDDENLP.upperBoxCons = [aDDENLP.upperBoxCons;certOptParamUB];
            else
                if aDDE.certOptParam.nVar > 0
                    
                    warning('upper boundaries for certain optimization parameters have to be handed over as a column vector with as many entries as parameters in the DDE. Using a vector of Inf instead')
                end
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
        
        
    % ======================================================================
    %> @brief intialize steady state constraint within a DDENLP instance
    %>
    %> @param aDDENLP instance of DDENLP
    % ====================================================================== 
        
        function initializeStSt( aDDENLP )
            % the initial guess for the steady state is corrected to find
            % an actual steady state
            options=aDDENLP.optionsInitEqCons;
            aDDENLP.stStCon = initStStConstraint(aDDENLP.stStCon,options);
        end
        
   	% ======================================================================
    %> @brief add a normal vector constraint to the optimization Problem
    %>
    %> @param aDDENLP instance of DDENLP
    %> @param type a string containing the requested manifold type
    %> @param augSynVSysHandlesHandle function handle describing the critical manifold
    %> @param nVSysHandle function handle describing the normal vectors of the critical manifold
    %> @param xGuess guess for critical state, instance of VariableVector
    %> @param alphaGuess guess for critical parameters, instance of VariableVector
    %> @param pGuess current certain parameters,  instance of VariableVector
    % ====================================================================== 
        
        function addNVCon( aDDENLP, type, augSysHandle, nVSysHandle, xGuess, alphaGuess, pGuess )
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
        
    % ======================================================================
    %> @brief initializes normal vector constraints by using a series of
    %>        methods of the class NVConstraint
    %>
    %> @param aDDENLP instance of DDENLP
    % ====================================================================== 
        
        function initNVCons( aDDENLP )
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
        
        
    % ======================================================================
    %> @brief moves the nominal point away from known critical manifolds
    %>        methods of the class NVConstraint
    %>
    %> @param aDDENLP instance of DDENLP
    %> @param steppingFactor (optional) step width of moving nominal point.
    %>        Smaller value means slower convergence but better numerics, default is 0.7
    %> @param distanceFactor (optional) desired final minimal distance
    % >       compared to property minDist. Default is 1.8
    %> @param iterations (optional) number of iterations. Default is 20
    % ====================================================================== 
        
        function moveAwayFromManifolds( aDDENLP, steppingFactor, distanceFactor,...
                iterations, varargin )
            % shifts nominal point further away from critical manifolds
            % this process works iteratively
     
            
            % if not handed over, assign default values
            if nargin==1
                steppingFactor = 0.7;
                distanceFactor = 1.8;
                iterations = 20;
            else
                if distanceFactor<1
                    warning('distanceFactor > 1 recommended ')
                end
            end
            
            refernceDistance = distanceFactor*aDDENLP.minDist*sqrt(aDDENLP.nAlpha);
            
            distance = NaN(1,length(aDDENLP.NVCon));            
            normalVectors = NaN(aDDENLP.nAlpha,length(aDDENLP.NVCon));
            
            
            for jj = 1:length(aDDENLP.NVCon)
                distance(1,jj) = aDDENLP.vars.critical(jj).l.values;
                normalVectors(:,jj) = aDDENLP.vars.critical(jj).r.values;
%                 quiver(aDDENLP.vars.critical(jj).alpha.values(1),aDDENLP.vars.critical(jj).alpha.values(2),...
%                 normalVectors(1,jj),normalVectors(2,jj),'k')
            end
            
            
            
            for ii=1:iterations
                % leave for loop when far away from critical manifolds
                if all(distance>=aDDENLP.minDist*sqrt(aDDENLP.nAlpha))
%                     plot(newNomPoint(1),newNomPoint(2),'k.');                    
%                     [uncertRegionAlpha, uncertRegionGamma] = circle( sqrt(2), aDDENLP.vars.nominal.alpha.values(1), aDDENLP.vars.nominal.alpha.values(2) );
%                     plot(uncertRegionAlpha, uncertRegionGamma);
                    fprintf('\nstopping to move the nominal point, %i iterations necessary to move nominal point sufficiently far away from critical maninfolds\n',ii-1);
                    break
                end
                newNomPoint = aDDENLP.vars.nominal.alpha.values +...
                    steppingFactor*normalVectors*...
                    max(zeros(1,length(aDDENLP.NVCon)),(refernceDistance-distance))';
%                 plot(newNomPoint(1),newNomPoint(2),'b.');
%                 drawnow        
                aDDENLP.vars.nominal.alpha.values = newNomPoint;
                
                %% update 
                aDDENLP.initializeStSt();
                aDDENLP.initNVCons();
                for jj = 1:length(aDDENLP.NVCon)
                    distance(1,jj) = aDDENLP.vars.critical(jj).l.values;
                    normalVectors(:,jj) = aDDENLP.vars.critical(jj).r.values;
                end
            end
        end
        
        
    % ======================================================================
    %> @brief evaluates the status of DDENLP
    %>
    %> @param aDDENLP instance of DDENLP
    % ====================================================================== 
        
        function evaluateStatus( aDDENLP )
            % checks the status of all existing constraints and assign the
            % status of the whole optimization problem
            if aDDENLP.stStCon.status==1
                aDDENLP.status=min(arrayfun(@(x)x.status,aDDENLP.NVCon));
            end
        end
        
        
        
        
    % ======================================================================
    %> @brief concatenates the various constraints in this DDENLP instance
    %>        for later optimitzation
    %>
    %> @param aDDENLP instance of DDENLP
    % ====================================================================== 
        function concatConstraints( aDDENLP )
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
        
        
        
    % ======================================================================
    %> @brief compares connection and normal vector
    %>
    %> @param aDDENLP instance of DDENLP
    %> @param indexCrit index of critical point for which comparison takes place
    % ====================================================================== 
        function compareConnectionNV(aDDENLP, indexCrit )
            % compares connection and normal vector
            
            connection = aDDENLP.vars.nominal.alpha.values - aDDENLP.vars.critical(indexCrit).alpha.values;
            connection = connection./norm(connection,2);
            
            r = aDDENLP.vars.critical.r.values;
            r = r/norm(r,2);
            disp([connection,r]);
                      
        end
        
        
    % ======================================================================
    %> @brief concatenates the various variables in this DDENLP instance
    %>        for later optimitzation
    %>
    %> @param aDDENLP instance of DDENLP
    % ====================================================================== 
        function concatInitPoints( aDDENLP )
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
        
    % ======================================================================
    %> @brief check if constraints hold for inital point
    %>
    %> @param aDDENLP instance of DDENLP
    % ====================================================================== 
        function checkConstraints( aDDENLP )
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
        
        
    	% ======================================================================
    %> @brief run optimization of DDENLP instance
    %>
    %> @param aDDENLP instance of DDENLP
    %> @param userDefinedOptions (optional) user defined optimization
    %>        options    
    %> @param Aineq (optional) user defined linear inequality constraints
    %> @param bineq (optional) user defined linear inequality constraints
    % ====================================================================== 
        
        function runOptim( aDDENLP, userDefinedOptions, Aineq, bineq, varargin )
            % this function starts the main optimization without
            % intermediate stability checks
            
            if nargin > 1
                options = userDefinedOptions;
            else
                options = aDDENLP.optionsMainOptim;
            end
            
            if nargin<=2
                Aineq=[];
                bineq=[];
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
                    fprintf('\n values of  inequality constraint values (which should be < 0)are: \n %f \n',c);
                    warning('nonlinear inequality constraints violated at initial point')
                end
            end
            
            [maxCeq,position] = max(abs(ceq));
            if (maxCeq>options.TolCon) && (~strcmp(callerFunction(2).name, 'DDENLP.runOptimWithStabChecks'))
                fprintf('%d\n',ceq);
                warning('nonlinear equality constraints violated at initial point for equality constraint no. %i, i.e. %d (measured by TolCon = %d).',position,abs(ceq(position)),(options.TolCon))
            end
            
            % run optimization
            
            [aDDENLP.optimVal,aDDENLP.optJ,aDDENLP.exitflag,aDDENLP.optimOutput,aDDENLP.lambda] = ...
                fmincon(J,x0,Aineq,bineq,Aeq,beq,lb,ub,nonlinCon,options);
            
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
        
    % ======================================================================
    %> @brief run optimization of DDENLP instance with intermediate
    %>        stability checks. Less iterations between stability checks
    %>        lead to slower optimization, because optimizer has to cold
    %>        start more frequently
    %>
    %> @param aDDENLP instance of DDENLP
    %> @param nIterBetweenStabChecks number of iterations between stability
    %>        checks
    %>    
    %> @return init initial point of  optimization
    %> @return final final point of optimization 
    %> @return maxEig biggest eigenvalue real part after optimization
    %> @return eigs rightmost eigenvalues after stability loss
    % ====================================================================== 
        
        function [init, final, maxEig, eigs] = runOptimWithStabChecks( aDDENLP, nIterBetweenStabChecks )
            % this function starts the main optimization with
            % intermediate stability checks
            
            options = aDDENLP.optionsMainOptim;
            maxIter = options.MaxIter;
            
            options.MaxIter = nIterBetweenStabChecks;
            options.Display = 'off';%iter-detailed';
            
            for ii = 1:ceil(maxIter/nIterBetweenStabChecks)
                aDDENLP.concatInitPoints();
                aDDENLP.runOptim(options);
                init = aDDENLP.deconstructInit();
                final = aDDENLP.deconstructOptimum();
                cost = aDDENLP.optJ;
                %                 [maxEig,eigs] = aDDENLP.checkStabilityPoint('nominal');
                [~,eigs] = aDDENLP.checkStabilityPoint(aDDENLP.vars.nominal);
                maxEig=eigs(aDDENLP.allowedEigsInClosedRightHP+1);
                fprintf('%d iterations completed, current cost is %d, eigenvalue with biggest real part was %d\n',ii*nIterBetweenStabChecks,cost,maxEig)
                
                callerFunction = dbstack;
                if real(maxEig) > aDDENLP.maxAllowedRealPart
                    if ~strcmp(callerFunction(2).name, 'DDENLP.runOptimAddingNewManifolds')
                        warning('system lost desired stability properties during optimization. You can use the properies DDENLP.initVal and DDENLP.optimVal to find out were stability was lost')
                    else
                        fprintf('system lost desired stability properties during optimization. Adding a new critical manifold to optimization problem\n');
                    end
                    break
                end
                % leave loop if optimization finished
                if aDDENLP.exitflag > 0
                    aDDENLP.runOptim();
                    break
                end
                
            end
        end
        
        
                
    % ======================================================================
    %> @brief run optimization of DDENLP instance with intermediate
    %>        stability checks. Algorithm tries to add new critical
    %>        manifold if one was crossed
    %>
    %> @param aDDENLP instance of DDENLP
    %> @param nIterBetweenStabChecks number of iterations between stability
    %>        checks and potential manifold adding
    % ====================================================================== 
        
        function runOptimAddingNewManifolds( aDDENLP, nIterBetweenStabChecks )
            
            for ii = 1:10
                [init,final,maxEig,~] = aDDENLP.runOptimWithStabChecks(nIterBetweenStabChecks);
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
                    
                    [handleManifold,handleNV] = aDDENLP.problemDDE.getHandles(type);
                    intermediatePoint = aDDENLP.findManifoldPointOnLine(type, handleManifold, init, final);
                    
                    %add new NVCON
                    aDDENLP.addNVCon(type,handleManifold,handleNV,intermediatePoint.x,intermediatePoint.alpha,init.p)
                    fprintf('\nadded new %s manifold', type);
                    
                    % and initialize it:
                    aDDENLP.deconstructInit('nom');
                    aDDENLP.concatInitPoints;                    
                    
                    aDDENLP.NVCon(end).findManifoldPoint(aDDENLP.NVCon(end).vars);
                    aDDENLP.NVCon(end).findClosestCriticalPoint(aDDENLP.vars.nominal.alpha);
                    aDDENLP.NVCon(end).findNormalVector(aDDENLP.vars.nominal.alpha);
                    aDDENLP.NVCon(end).findConnection(aDDENLP.vars.nominal.alpha);
                    
                    aDDENLP.concatConstraints;
                else
                    break
                end
            end
        end
        
                        
    % ======================================================================
    %> @brief look for a manifold point between two given points
    %>
    %> @param ~ (ignore first input)
    %> @param type string of expected manifold type
    %> @param manifoldHandle function handle of manifold
    %> @param point1 first point
    %> @param point2 second point
    %> 
    %> @return intermediate Point
    % ====================================================================== 
        
        
        
        function  intermediatePoint = findManifoldPointOnLine( ~, type, manifoldHandle, point1, point2 )
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
            
            options=optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',2000,'MaxFunEvals',2000000,'display','off','TolFun',1e-7,'TolX',1e-9);
            
            [y,res,localexitflag] = fsolve(rhs,y0,options);
            if localexitflag > 0
                fprintf('\nlooked for a %s point on the line between optimization points and appearantly found one, exitflag was %d\n', type, localexitflag)
            else
                disp(res);
                warning('did not find critical point (%s) on the given line, fsolve exitflag was %d',type,localexitflag)
            end

            x = y(1:point1.x.nVar);
            mu = y(point1.x.nVar+1);
            if (mu < 0) || (mu > 1)
                warning('converged to a point that is not between the given points, mu = %d', mu)
            end
            alpha = alphaInterp( mu );
            p =  pInterp( mu );
            omega = y(point1.x.nVar+2);
            w1 = y(point1.x.nVar+3:2*point1.x.nVar+2);
            w2 = y(2*point1.x.nVar+3:3*point1.x.nVar+2);
            
            intermediatePoint.x = VariableVector(x,Inf,point1.x.names);
            intermediatePoint.alpha = VariableVector(alpha, Inf,point1.alpha.names);
            intermediatePoint.p = VariableVector(p, Inf,point1.p.names);
            intermediatePoint.omega = VariableVector(omega, Inf);
            intermediatePoint.w1 = VariableVector(w1, Inf);
            intermediatePoint.w2 = VariableVector(w2, Inf);
        end
        
        
        %% FROM HERE UPWARDS: ENDING MARK
        
                                
    % ======================================================================
    %> @brief reconstruct inital variables from single vector
    %>
    %> @param aDDENLP instance of DDENLP class
    %> @param type (optional) string describing which points to reconstruct
    %> 
    %> @return (optional) return reconstruction to external variable
    % ====================================================================== 
        
        function [varargout] = deconstructInit(aDDENLP, type, varargin)
            
            init = aDDENLP.initVal;
            outputVars = aDDENLP.vars;
            
            if nargin == 1
                type='all';
            end
            
            if strcmp(type,'all') || strcmp(type,'nom')
                outputVars.nominal.x.values = init(outputVars.nominal.x.index);
                outputVars.nominal.alpha.values = init(outputVars.nominal.alpha.index);
                outputVars.nominal.p.values = init(outputVars.nominal.p.index);
            end
            
            
            if strcmp(type,'all') || strcmp(type,'crit')
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
            end
            
            if nargout > 0
                varargout{1} = varCollection( 'nominal', NaN, outputVars.nominal.x, outputVars.nominal.alpha, outputVars.nominal.p);
            else
                aDDENLP.vars = outputVars;
            end
        end
        
    % ======================================================================
    %> @brief reconstruct optimal variables from single vector
    %>
    %> @param aDDENLP instance of DDENLP class
    %> 
    %> @return (optional) return reconstruction to external variable
    % ====================================================================== 
        
        function varargout = deconstructOptimum(aDDENLP)
            
            
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
            if nargout > 0
                varargout{1} =  varCollection( 'nominal', NaN, outputVars.nominal.x, outputVars.nominal.alpha, outputVars.nominal.p);
            else
                aDDENLP.vars = outputVars;
            end
            
        end
        
    % ======================================================================
    %> @brief calculate eigenvalues at a given point in parameter space
    %>
    %> @param aDDENLP instance of DDENLP class
    %> @param point string desribing at which point type the stability shall be evaluated
    %> 
    %> @return maxRealPart maximal real part of all eigenvalues
    %> @return eigs some rightmost eigenvalues
    % ====================================================================== 
        
        function [ maxRealPart, eigs ] = checkStabilityPoint( aDDENLP, point )
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
            
            for ii = 1 : size(point,1)
                [maxRealPartTMP, xCorr, eigs] = checkStability(funcs,point(ii).alpha.values,point(ii).x.values,aDDENLP.numMinEig);
                if isempty(maxRealPartTMP)
                    error('did not find an eigevalue with Re(lambda)>aDDENLP.numMinEig. Please chose a smaller value for aDDENLP.numMinEig')
                end
                maxRealPart(ii) = maxRealPartTMP;
                if exist('type','var')
                    switch type
                        case 'nominal'
                            aDDENLP.vars.nominal(ii,1).maxEig = maxRealPart(ii);
                            aDDENLP.vars.nominal(ii,1).eigs = eigs;
                        case 'critical'
                            aDDENLP.vars.critical(ii,1).maxEig = maxRealPart(ii);
                            aDDENLP.vars.critical(ii,1).eigs = eigs;
                        case 'vertex'
                            aDDENLP.vars.vertex(ii,1).x.values = xCorr;
                            aDDENLP.vars.vertex(ii,1).maxEig = maxRealPart(ii);
                            aDDENLP.vars.vertex(ii,1).eigs = eigs;
                        case 'random'
                            aDDENLP.vars.random(ii,1).x.values = xCorr;
                            aDDENLP.vars.random(ii,1).maxEig = maxRealPart(ii);
                            aDDENLP.vars.random(ii,1).eigs = eigs;
                    end
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
        
    % ======================================================================
    %> @brief calculate stability an vertices of uncertianty region
    %>
    %> @param aDDENLP instance of DDENLP class
    %> 
    %> @return maxRealPart maximal real part of all eigenvalues
    % ====================================================================== 
        
        function maxRealPart = checkStabilityAtVertices( aDDENLP )
            
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
        
    % ======================================================================
    %> @brief calculate stability at random points within uncertianty region
    %>
    %> @param aDDENLP instance of DDENLP class
    %> @param baseForNumberOfPoints base for the number of random points
    %> @param seedForRandomNumbers (optional) seed for pseudo random number
    %>        generator
    %>
    %> @return maxRealPart maximal real part of all eigenvalues
    % ====================================================================== 
        
        function maxRealPart = checkStabilityAtRandom( aDDENLP, baseForNumberOfPoints, seedForRandomNumbers, varargin )
            % calculates max(real(eigenvalues)) for random points within 
            % the uncertainty region
            
            if nargin > 2
                rng(seedForRandomNumbers);
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
        
        
    % ======================================================================
    %> @brief run a simulation of the optimum. Overloading of ddesd
    %>
    %> @param aDDENLP instance of DDENLP class
    %> @param point parameters to use for simulation
    %> @param history for simulation (like inital point ODE, but for DDE)
    %> @param tspan time spaned by simulation
    %> @param options options for solver
    %>
    %> @return sol is solution struct as known from ode45 etc
    % ====================================================================== 
        
        function sol = ddesd( aDDENLP, point, history, tspan, options )
            
            myDDE=@(~,x,xtau)aDDENLP.problemDDE.rhs(x,xtau,point.alpha.values,point.p.values)';
            myDelays=@(t,x)t-aDDENLP.problemDDE.delays(x,point.alpha.values,point.p.values)';
            
            sol = ddesd(myDDE,myDelays,history,tspan,options) ;
        end
        
        %         export2DDEBIFTOOL % nicht bei class DDE, damit Box
        %         constraints auch übergeben werden können
        
    end
    
end

