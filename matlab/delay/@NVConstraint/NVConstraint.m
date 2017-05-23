classdef NVConstraint < EqualityConstraint
    %NVCONSTRAINT the instances of this class are objects representing
    %normal vector constraints. The procedures initialize those constraints
    
    properties
        
        type % string
        
        eqAugSys % EqualityConstraint object
        eqNVSys % EqualityConstraint object
        eqConnect % EqualityConstraint object
        
        nVarAugSys % number of variables in augemented system
        nVarNVSys % number of variables in normal vector system
        
        inequalities % function handle
        optionsEqConsInit % options for numerical solver
        optionsInitOptim % options for auxiliary optimization
        
    end
    
    
    properties(SetAccess=protected)
        inequalityIndex % index of the inequality constraint of this NVConstraint within all constaints
        problemDDE % differential equation for calculation of stability
        numMinEig % numerical constant for the ODE approximation of the DDE
    end
    
    methods
        function aNVCon = NVConstraint(aDDENLP,type,augSysHandle,nVSysHandle,nVvars) % constructor
            
            %% translate human readable input format to solver readable format
            switch type
                case 'fold'
                    nEqAugSys=2*aDDENLP.vars.nominal.x.nVar+1;
                    nEqNVSys=aDDENLP.vars.nominal.alpha.nVar;
                    augSysHandle=@(y)augSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.w1.index))';
                    nVSysHandle=@(y)nVSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.w1.index),y(nVvars.r.index))';
                case 'modfold'
                    nEqAugSys=2*aDDENLP.vars.nominal.x.nVar+1;
                    nEqNVSys=2*aDDENLP.vars.nominal.x.nVar+aDDENLP.vars.nominal.alpha.nVar+1;
                    augSysHandle=@(y)augSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.w1.index))';
                    nVSysHandle=@(y)nVSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.w1.index),y(nVvars.v1.index),y(nVvars.g1.index),y(nVvars.u.index),y(nVvars.r.index))';
                case 'hopf'
                    nEqAugSys=3*aDDENLP.vars.nominal.x.nVar+2;
                    nEqNVSys=3*aDDENLP.vars.nominal.x.nVar+aDDENLP.vars.nominal.alpha.nVar+2;
                    augSysHandle=@(y)augSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.omega.index),y(nVvars.w1.index),y(nVvars.w2.index))';
                    nVSysHandle=@(y)nVSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.omega.index),y(nVvars.w1.index),y(nVvars.w2.index),y(nVvars.v1.index),y(nVvars.v2.index),y(nVvars.g1.index),y(nVvars.g2.index),y(nVvars.u.index),y(nVvars.r.index))';
                case 'modhopf'
                    nEqAugSys=3*aDDENLP.vars.nominal.x.nVar+2;
                    nEqNVSys=3*aDDENLP.vars.nominal.x.nVar+aDDENLP.vars.nominal.alpha.nVar+2;
                    augSysHandle=@(y)augSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.omega.index),y(nVvars.w1.index),y(nVvars.w2.index))';
                    nVSysHandle=@(y)nVSysHandle(y(nVvars.x.index),y(nVvars.alpha.index),y(nVvars.p.index),y(nVvars.omega.index),y(nVvars.w1.index),y(nVvars.w2.index),y(nVvars.v1.index),y(nVvars.v2.index),y(nVvars.g1.index),y(nVvars.g2.index),y(nVvars.u.index),y(nVvars.r.index))';
                otherwise
                    error(['invalid manifold type ',type,', could not construct normal vector constraint object'])
            end
            
            %% define connection Constraint
            connectionConstraint=@(y)...
                y(aDDENLP.vars.nominal.alpha.index)...
                -y(nVvars.alpha.index)...
                -y(nVvars.l.index)*y(nVvars.r.index)/norm(y(nVvars.r.index),2);
            
            %% write constraints in class properties
            aNVCon=aNVCon@EqualityConstraint(@(y)[augSysHandle(y);nVSysHandle(y);connectionConstraint(y)],...
                nEqAugSys+nEqNVSys+aDDENLP.vars.nominal.alpha.nVar,...
                nVvars,...
                aDDENLP.occupiedEqs);

            %% define solver options
            aNVCon.optionsEqConsInit= aDDENLP.optionsInitEqCons;
            aNVCon.optionsInitOptim = aDDENLP.optionsInitOptim;
            
            aNVCon.type=type;
            aNVCon.eqAugSys=augSysHandle;
            aNVCon.eqNVSys=nVSysHandle;
            
            
            aNVCon.eqAugSys=EqualityConstraint(augSysHandle,nEqAugSys,nVvars,aDDENLP.occupiedEqs);
            aNVCon.eqNVSys=EqualityConstraint(nVSysHandle,nEqNVSys,nVvars,aDDENLP.occupiedEqs+nEqAugSys);
            aNVCon.eqConnect=EqualityConstraint(connectionConstraint,nVvars.alpha.nVar,nVvars.nVar,aDDENLP.occupiedEqs+nEqAugSys+nEqNVSys);
            %             aNVCon.eqConnect=EqualityConstraint(connectionConstraint,1,nVvars.nVar,aDDENLP.occupiedEqs+nEqAugSys+nEqNVSys);
            
            
            aNVCon.inequalities = @(y)aDDENLP.minDist*sqrt(aDDENLP.nAlpha)-y(aNVCon.vars.l.index);
            
            aNVCon.nVarAugSys=nVvars.x.nVar+nVvars.alpha.nVar+nVvars.omega.nVar+nVvars.w1.nVar+nVvars.w2.nVar;
            aNVCon.nVarNVSys=nVvars.v1.nVar+nVvars.v2.nVar+nVvars.g1.nVar+nVvars.g2.nVar+nVvars.u.nVar+nVvars.r.nVar;
            
            aNVCon.inequalityIndex=aDDENLP.occupiedIneqs+1;
            
            aNVCon.problemDDE = aDDENLP.problemDDE;
            aNVCon.numMinEig = aDDENLP.numMinEig;
            
        end
        
        function prepareInitialGuess(aNVCon,aVarCollection)
            % rotates complex eigenvector to make real and imaginary part
            % orthogonal
            % cf. Proof of Lemma 1 in
            % [https://doi.org/10.1109/CDC.2016.7798469]
            
            % get real and imaginary part
            realW = aVarCollection.w1.values;
            imagW = aVarCollection.w2.values;
            
            w = realW +1i*imagW;            
            phi = 0;
            
            if ~isreal(w)
                % find rotational angle
                f=@(PHI)cos(2*PHI)*real(w)'*imag(w) + 0.5*sin(2*PHI)*(real(w)'*real(w)-imag(w)'*imag(w)); 
                phi=fsolve(f,phi,aNVCon.optionsEqConsInit);
            end
            % make the transform
            w = exp(1i*phi)*w;
            
            % reassign numerical values
            aVarCollection.w1.values = real(w);
            aVarCollection.w2.values = imag(w);
            
        end
        
        function findManifoldPoint(aNVCon,aVarCollection)
            % tries to find a point on the critical manifold
            
            % evaluate number of leading variables, that are not part of this
            % constraint
            offset=aNVCon.vars.x.index(1)-1;
            
            otherVariables=zeros(offset,1);
            otherVariables(aVarCollection.p.index)=aVarCollection.p.values;
            
            rhs = @(y)aNVCon.eqAugSys.conFun([otherVariables;y]);
            
            if nargin>1
                x0=[aVarCollection.x.values;
                    aVarCollection.alpha.values;
                    aVarCollection.omega.values;
                    aVarCollection.w1.values;
                    aVarCollection.w2.values;];
            else
                fprintf('using default initial guess in manifold point search\n')
                x0 = ones(aNVCon.nVarAugSys,1);
            end
            
            [x,res,exitflag] = fsolve(rhs,x0,aNVCon.optionsEqConsInit);
            
            aNVCon.vars.x.values = x(aNVCon.vars.x.index-offset);
            aNVCon.vars.alpha.values = x(aNVCon.vars.alpha.index-offset);
            aNVCon.vars.omega.values = x(aNVCon.vars.omega.index-offset);
            aNVCon.vars.w1.values = x(aNVCon.vars.w1.index-offset);
            aNVCon.vars.w2.values = x(aNVCon.vars.w2.index-offset);
            
            if exitflag > 0
                aNVCon.status = 2;
                
                callerFunction = dbstack;
                if (length(callerFunction)>2) && ~strcmp(callerFunction(3).name, 'DDENLP.moveAwayFromManifolds')
                    fprintf('\nfound point on critical manifold, fsolve exitflag was %d\n', exitflag)
                end
            else
                warning('did not find critical point, fsolve exitflag was %d',exitflag)
                disp(res);
            end
        end
        
        function findClosestCriticalPoint(aNVCon,alphaNom)
            
            % check current status
            if aNVCon.status<2
                error('cannot search for closest critical point without critical point');
            end
            
            offset = aNVCon.vars.x.index(1)-1;
            
            otherVariables=zeros(offset,1);
            otherVariables(aNVCon.vars.p.index)=aNVCon.vars.p.values;
            
            manifoldCon = @(y)deal([],aNVCon.eqAugSys.conFun([otherVariables;y]));
            
            x0 = [aNVCon.vars.x.values;...
                aNVCon.vars.alpha.values;...
                aNVCon.vars.omega.values;...
                aNVCon.vars.w1.values;...
                aNVCon.vars.w2.values];
            
            % perform actual search for closest critical point
            J=@(y)norm(alphaNom.values-y(aNVCon.vars.alpha.index-offset),2);
            [x,l,exitflag] = fmincon(J,x0,[],[],[],[],[],[],manifoldCon,aNVCon.optionsInitOptim);
            
            % extract values
            aNVCon.vars.x.values = x(aNVCon.vars.x.index-offset);
            aNVCon.vars.alpha.values = x(aNVCon.vars.alpha.index-offset);
            aNVCon.vars.omega.values = x(aNVCon.vars.omega.index-offset);
            aNVCon.vars.w1.values = x(aNVCon.vars.w1.index-offset);
            aNVCon.vars.w2.values = x(aNVCon.vars.w2.index-offset);
            aNVCon.vars.l.values = l;
            
            % give feedback
            if exitflag>0
                aNVCon.status=3;
                callerFunction = dbstack;
                if (length(callerFunction)>2) && ~strcmp(callerFunction(3).name, 'DDENLP.moveAwayFromManifolds')
                    fprintf('found closest point on critical manifold at distance l = %f, fmincon exitflag was %d\n',l,exitflag)
                end
            else
                warning('did not find closest critical point, fmincon exitflag was %d', exitflag)
            end
        end
        
        function  findNormalVector(aNVCon,alphaNom)
            
            directionMode = 2;
            
            if aNVCon.status < 2
                error('cannot search for normal vector without a critical point');
            end
            
            offset=aNVCon.vars.x.index(1)-1;
            
            otherVariables=zeros(offset,1);
            otherVariables(aNVCon.vars.p.index)=aNVCon.vars.p.values;
            
            maniPoint = [aNVCon.vars.x.values;...
                aNVCon.vars.alpha.values;...
                aNVCon.vars.omega.values;...
                aNVCon.vars.w1.values;...
                aNVCon.vars.w2.values];
            
            rhs = @(y)aNVCon.eqNVSys.conFun([otherVariables; maniPoint; y]);
            
            x0 = ones(aNVCon.nVarNVSys,1);
            
            [x,res,exitflag] = fsolve(rhs,x0,aNVCon.optionsEqConsInit);
            
            if ~strcmp(aNVCon.type,'fold')
                offset=aNVCon.vars.v1.index(1)-1;
            else
                offset=aNVCon.vars.w1.index(end);
            end
            %% THIS IS A TEMPORARY CODE SEGMENT WHICH ASSUMES THE NOMINAL POINT IS IN THE STABLE REGION
            % IT FLIPS THE NORMAL VECTOR ALLWAYS TOWARDS THE NOMINAL POINT. LATER THE NORMAL VECTOR HAS TO POINT TO THE STABLE REGION
            
            switch directionMode
                case 1      % FLIPS THE NORMAL VECTOR ALLWAYS TOWARDS THE NOMINAL POINT
                    
                    r = x(aNVCon.vars.r.index-offset);
                    r = r/norm(r,2);
                    alphaNom = alphaNom.values;
                    alphaCrit = aNVCon.vars.alpha.values;
                    direction = sign((r'*r)\r'*(alphaNom-alphaCrit));
                    
                case 2  % FLIPS THE NORMAL VECTOR ALLWAYS TOWARDS THE STABLE REGION, determined by eigenvalues
                    
                    shiftLength = 1e-3;
                    
                    r = x(aNVCon.vars.r.index-offset);
                    r = r/norm(r,2);
                    
                    alpha1 = aNVCon.vars.alpha.values - shiftLength*r;                    
                    alpha2 = aNVCon.vars.alpha.values + shiftLength*r;
                    
                    xStSt = aNVCon.vars.x.values;
                    
                    %% prepare stability calculations
                    p=aNVCon.vars.p.values;
                    
                    myDDE = @(xx,alpha)aNVCon.problemDDE.rhs(xx(:,1),xx(:,2:end),alpha,p)';
                    myNtau = aNVCon.problemDDE.ntau;
                    myDelays = @(k,xx,alpha)(k==(1:myNtau))*aNVCon.problemDDE.delays(xx,alpha,p);
                    
                    funcs = set_funcs('sys_rhs',myDDE,...
                        'sys_ntau',myNtau,...
                        'sys_tau',myDelays);
                    
                    [lambda1,~,~] = checkStability(funcs,alpha1,xStSt,aNVCon.numMinEig);
                    if isempty(lambda1)
                        error('did not find an eigevalue with Re(lambda)>aDDENLP.numMinEig. Please chose a smaller value for aDDENLP.numMinEig')
                    end
                    
                    [lambda2,~,~] = checkStability(funcs,alpha2,xStSt,aNVCon.numMinEig);
                    if isempty(lambda2)
                        error('did not find an eigevalue with Re(lambda)>aDDENLP.numMinEig. Please chose a smaller value for aDDENLP.numMinEig')
                    end
                    
                    %% evaluate
                    direction = sign(lambda1-lambda2);
                    if direction == 0
                        warning('could not determine which region is stable and which is unstable')
                    end
                    %                 case 3  % FLIPS THE NORMAL VECTOR ALLWAYS TOWARDS THE STABLE REGION, determined by analytical expression
            end
            
            % END OF THE TEMPORARY SEGMENT
            
            
            aNVCon.vars.v1.values = direction*x(aNVCon.vars.v1.index-offset);
            aNVCon.vars.v2.values = direction*x(aNVCon.vars.v2.index-offset);
            aNVCon.vars.g1.values = direction*x(aNVCon.vars.g1.index-offset);
            aNVCon.vars.g2.values = direction*x(aNVCon.vars.g2.index-offset);
            aNVCon.vars.u.values = direction*x(aNVCon.vars.u.index-offset);
            aNVCon.vars.r.values = direction*x(aNVCon.vars.r.index-offset);
            
            if exitflag>0
                if aNVCon.status == 3
                    aNVCon.status=4;
                    callerFunction = dbstack;
                    if (length(callerFunction)>2) && ~strcmp(callerFunction(3).name, 'DDENLP.moveAwayFromManifolds')
                        fprintf('found normal vector on closest critical manifold point, fsolve exitflag was %d\n',exitflag)
                    end
                else
                    fprintf('found normal vector on a critical manifold point, fsolve exitflag was %d\n',exitflag)
                end
                
            else
                disp(res)
                warning('did not find normal vector on critical point, fsolve exitflag was %d\n', exitflag)
            end
        end
        
        function findConnection(aNVCon,alphaNom)
            
            if aNVCon.status<4
                error('cannot search for connection without normal vector');
            end
            
            offset=aNVCon.vars.x.index(1)-1;
            
            x0 = [aNVCon.vars.x.values;...
                aNVCon.vars.alpha.values;...
                aNVCon.vars.omega.values;...
                aNVCon.vars.w1.values;...
                aNVCon.vars.w2.values;...
                aNVCon.vars.v1.values;...
                aNVCon.vars.v2.values;...
                aNVCon.vars.g1.values;...
                aNVCon.vars.g2.values;...
                aNVCon.vars.u.values;...
                aNVCon.vars.r.values;...
                aNVCon.vars.l.values];
            
            otherVariables=zeros(offset,1);
            otherVariables(alphaNom.index)=alphaNom.values;
            otherVariables(aNVCon.vars.p.index)=aNVCon.vars.p.values;
            
            rhs = @(y)aNVCon.conFun([otherVariables;y]);
            
            [x,res,exitflag] = fsolve(rhs,x0,aNVCon.optionsEqConsInit);
            
            aNVCon.vars.x.values = x(aNVCon.vars.x.index-offset);
            aNVCon.vars.alpha.values = x(aNVCon.vars.alpha.index-offset);
            aNVCon.vars.omega.values = x(aNVCon.vars.omega.index-offset);
            aNVCon.vars.w1.values = x(aNVCon.vars.w1.index-offset);
            aNVCon.vars.w2.values = x(aNVCon.vars.w2.index-offset);
            aNVCon.vars.v1.values = x(aNVCon.vars.v1.index-offset);
            aNVCon.vars.v2.values = x(aNVCon.vars.v2.index-offset);
            aNVCon.vars.g1.values = x(aNVCon.vars.g1.index-offset);
            aNVCon.vars.g2.values = x(aNVCon.vars.g2.index-offset);
            aNVCon.vars.u.values = x(aNVCon.vars.u.index-offset);
            aNVCon.vars.r.values = x(aNVCon.vars.r.index-offset);
            aNVCon.vars.l.values = x(aNVCon.vars.l.index-offset);
            
            if exitflag>0
                aNVCon.status=5;
                callerFunction = dbstack;
                if ~strcmp(callerFunction(3).name, 'DDENLP.moveAwayFromManifolds')
                    fprintf('found connection of nominal and critical point, fsolve exitflag was %d \n\n',exitflag)
                end
            else
                disp(res)
                warning('no connection found, fsolve exitflag was %d',exitflag)
            end
        end
    end
    
end

