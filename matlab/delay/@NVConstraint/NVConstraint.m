classdef NVConstraint < EqualityConstraint
    %NVCONSTRAINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        type % string
        
        eqAugSys % EqualityConstraint object
        eqNVSys % EqualityConstraint object
        eqConnect % EqualityConstraint object
        
        nVarAugSys
        nVarNVSys
        
        inequalities % function handle
        optionsEqConsInit
        optionsInitOptim
        
    end
    
    
    properties(SetAccess=protected)
        inequalityIndex
    end
    
    methods
        function aNVCon = NVConstraint(aDDENLP,type,augSysHandle,nVSysHandle,nVvars)
           
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
             
            
            aNVCon.inequalities = @(y)aDDENLP.minDist*sqrt(aDDENLP.nAlpha)-y(aNVCon.vars.l.index);
            
            aNVCon.nVarAugSys=nVvars.x.nVar+nVvars.alpha.nVar+nVvars.omega.nVar+nVvars.w1.nVar+nVvars.w2.nVar;
            aNVCon.nVarNVSys=nVvars.v1.nVar+nVvars.v2.nVar+nVvars.g1.nVar+nVvars.g2.nVar+nVvars.u.nVar+nVvars.r.nVar;
                       
            aNVCon.inequalityIndex=aDDENLP.occupiedIneqs+1;
            
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
            
            [x,~,exitflag] = fsolve(rhs,x0,aNVCon.optionsEqConsInit);
                        
            aNVCon.vars.x.values = x(aNVCon.vars.x.index-offset);
            aNVCon.vars.alpha.values = x(aNVCon.vars.alpha.index-offset);
            aNVCon.vars.omega.values = x(aNVCon.vars.omega.index-offset);
            aNVCon.vars.w1.values = x(aNVCon.vars.w1.index-offset);
            aNVCon.vars.w2.values = x(aNVCon.vars.w2.index-offset);
            
            if exitflag > 0
                aNVCon.status = 2;
                fprintf('\nfound point on critical manifold, fsolve exitflag was %d\n', exitflag)
            else
                warning('did not find critical point, fsolve exitflag was %d',exitflag)
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
                fprintf('found closest point on critical manifold at distance l = %f, fsolve exitflag was %d\n',l,exitflag)
              %  disp(aNVCon.vars.alpha.values')
            else
                 warning('did not find closest critical point, fmincon exitflag was %d', exitflag)
            end
        end
        
        function  findNormalVector(aNVCon,alphaNom)
            
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
            
            [x,~,exitflag] = fsolve(rhs,x0,aNVCon.optionsEqConsInit);
            
            offset=aNVCon.vars.v1.index(1)-1;
            
            %% THIS IS A TEMPORARY CODE SEGMENT WHICH ASSUMES THE NOMINAL POINT IS IN THE STABLE REGION
            % IT FLIPS THE NORMAL VECTOR ALLWAYS TOWARDS THE NOMINAL POINT. LATER THE NORMAL VECTOR HAS TO POINT TO THE STABLE REGION
            
            r=x(aNVCon.vars.r.index-offset);
            r=r/norm(r,2);
            alphaNom=alphaNom.values;
            alphaCrit=aNVCon.vars.alpha.values;
            direction=sign((r'*r)\r'*(alphaNom-alphaCrit));
            
            
            % END OF THE TEMPORARY SEGMENT

            aNVCon.vars.v1.values = direction*x(aNVCon.vars.v1.index-offset);
            aNVCon.vars.v2.values = direction*x(aNVCon.vars.v2.index-offset);
            aNVCon.vars.g1.values = direction*x(aNVCon.vars.g1.index-offset);
            aNVCon.vars.g2.values = direction*x(aNVCon.vars.g2.index-offset);
            aNVCon.vars.u.values = direction*x(aNVCon.vars.u.index-offset);
            aNVCon.vars.r.values = direction*x(aNVCon.vars.r.index-offset);
            
            
            %             blaNV=rhs([aNVCon.vars.v1.values;...
            %                 aNVCon.vars.v2.values;...
            %                 aNVCon.vars.g1.values;...
            %                 aNVCon.vars.g2.values;...
            %                 aNVCon.vars.u.values;...
            %                 aNVCon.vars.r.values;...
            %                 aNVCon.vars.l.values])
            %
            
            if exitflag>0
                if aNVCon.status == 3
                    aNVCon.status=4;
                    fprintf('found normal vector on closest critical manifold point, fsolve exitflag was %d\n',exitflag)
                else
                    fprintf('found normal vector on a critical manifold point, fsolve exitflag was %d\n',exitflag)
                end
               
            else
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
% 
%             l=aNVCon.vars.l.values
%             
%             bla=rhs(x0);
%             connectionResBefore=bla(end-11:end)
%             
%             r=aNVCon.vars.r.values;
%             r=r/norm(r,2);
%             l=aNVCon.vars.l.values;
%             alphaNom=alphaNom.values;
%             alphaCrit=aNVCon.vars.alpha.values;
%             
%             checkBefore=[l*r,alphaNom-alphaCrit]
%             
%           [bla(end-11:end)/norm(bla(end-11:end),2),aNVCon.vars.r.values/norm(aNVCon.vars.r.values,2)]

            
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
                fprintf('found connection of nominal and critical point, fsolve exitflag was %d \n',exitflag)
            else
                disp(res)
                warning('no connection found, fsolve exitflag was %d',exitflag)
            end
        end
    end
    
end

