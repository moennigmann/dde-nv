%> @brief The instances of this class generate visualization data
%> @author Jonas Otten
%> @date 18 Jul 2017

%>======================================================================
%> @brief The instances of this class generate visualization data
%
%>  this class is allows easy visualization of optimization results by implementing a 
%>  quasi-arclength numerical continuation of critical manifolds.
%>
%> @author Jonas Otten
%> @date 18 Jul 2017 ======================================================================



classdef ManifoldSlice < handle
    %this class is creates slices of critical manifolds for later plotting
    %   Detailed explanation goes here
    
    properties
        freeParamIndices
        nManiPoints = 50;
        initStepLength = 0.01;
        stepLength   % changes during excution
        point
        eqAugSys
        eqNVSys
        maxStepLength = 0.5
        showStepsFlag = 0
        
        lowerBoxCons
        upperBoxCons
        
        debugFlag = 0;
        
    end

    
    % ======================================================================
    %> @brief Class constructor
    %>
    %> This function constructs instances of the class ManifoldSlice
    %> @param aNVCon an instance of the class NVConstraint. Its critical manifold will be visualized
    %> @param continParamsInd the index of the parameters used for continuation
    %> @param varargin enables to call constructor ommitting continParamsInd
    %>
    %> @return instance of the ManifoldSlice class.
    % ======================================================================    
    
    
    methods
        function aManifoldSlice = ManifoldSlice(aNVCon, continParamsInd, varargin) % constructor
            % CONSTRUCTOR OF THE CLASS MANIFOLDSLICE
            
            if nargin < 2
                [~,orderingNVIndexes] = sort(abs(aNVCon.vars.r.values),'descend'); % find largest entries of normal vector
                aManifoldSlice.freeParamIndices = aNVCon.vars.r.index(orderingNVIndexes(1:2))-aNVCon.vars.r.index(1)+1;
            else
                aManifoldSlice.freeParamIndices = continParamsInd;
            end
            
            aManifoldSlice.eqAugSys = aNVCon.eqAugSys;
            aManifoldSlice.eqNVSys = aNVCon.eqNVSys;
            aManifoldSlice.stepLength = aManifoldSlice.initStepLength;
            
            point.x=aNVCon.vars.x.copy;
            point.alpha=aNVCon.vars.alpha.copy;
            point.p=aNVCon.vars.p.copy;
            point.omega=aNVCon.vars.omega.copy;
            point.w1=aNVCon.vars.w1.copy;
            point.w2=aNVCon.vars.w2.copy;
            point.v1=aNVCon.vars.v1.copy;
            point.v2=aNVCon.vars.v2.copy;
            point.g1=aNVCon.vars.g1.copy;
            point.g2=aNVCon.vars.g2.copy;
            point.u=aNVCon.vars.u.copy;
            point.r=aNVCon.vars.r.copy;
            point.l=aNVCon.vars.l.copy;
            
            point.exitflagNV = -Inf;
            point.exitflagManifold = -Inf;
            
            aManifoldSlice.point  = point;
            
            aManifoldSlice.lowerBoxCons = -Inf(point.alpha.nVar,1);
            aManifoldSlice.upperBoxCons = +Inf(point.alpha.nVar,1);
            
            
        end
        
        
    % ======================================================================
    %> @brief runs numerical continuation in both directions
    %>
    %> @param aManifoldSlice instance of this class
    %> @param n number of steps to take during continuation
    %>
    % ======================================================================
        
        function maniContin2DbothDirections ( aManifoldSlice, n )
            for jj = 1:length(aManifoldSlice)
                aManifoldSlice(jj).nManiPoints = n;
                aManifoldSlice(jj).manifoldContinuation2D( 1);
                aManifoldSlice(jj).point = aManifoldSlice(jj).point(end:-1:1);
                aManifoldSlice(jj).stepLength = aManifoldSlice(jj).initStepLength;
                aManifoldSlice(jj).manifoldContinuation2D(-1);
            end
        end
        
        
        
    % ======================================================================
    %> @brief runs numerical continuation in one directions
    %>
    %> @param aManifoldSlice instance of this class
    %> @param direction determines direction for continuation
    %>
    % ======================================================================
        function manifoldContinuation2D( aManifoldSlice, direction )
            
            for jj = length(aManifoldSlice)
                for ii = length(aManifoldSlice(jj).point) :length(aManifoldSlice(jj).point) + aManifoldSlice(jj).nManiPoints-1
                    
                    
                    %% find tangent space and normal space
                    tngntSpc = zeros(aManifoldSlice(jj).point(ii).r.nVar,1);
                    tngntSpc(aManifoldSlice(jj).freeParamIndices(1)) =...
                        + aManifoldSlice(jj).point(ii).r.values(aManifoldSlice(jj).freeParamIndices(2));
                    tngntSpc(aManifoldSlice(jj).freeParamIndices(2)) =...
                        - aManifoldSlice(jj).point(ii).r.values(aManifoldSlice(jj).freeParamIndices(1));
                    
                    tngntSpc = sign(direction)*tngntSpc/norm(tngntSpc,2);
                    
                    nrmlSpc = zeros(aManifoldSlice(jj).point(ii).r.nVar,1);
                    nrmlSpc(aManifoldSlice(jj).freeParamIndices(1)) = aManifoldSlice(jj).point(ii).r.values(aManifoldSlice(jj).freeParamIndices(1));
                    nrmlSpc(aManifoldSlice(jj).freeParamIndices(2)) = aManifoldSlice(jj).point(ii).r.values(aManifoldSlice(jj).freeParamIndices(2));
                    
                    nrmlSpc=nrmlSpc/norm(nrmlSpc,2);
                    
                    %% extract data of temporay point
                    tmpPoint.x=aManifoldSlice(jj).point(ii).x.copy;
                    tmpPoint.alpha=aManifoldSlice(jj).point(ii).alpha.copy;
                    tmpPoint.p=aManifoldSlice(jj).point(ii).p.copy;
                    tmpPoint.omega=aManifoldSlice(jj).point(ii).omega.copy;
                    tmpPoint.w1=aManifoldSlice(jj).point(ii).w1.copy;
                    tmpPoint.w2=aManifoldSlice(jj).point(ii).w2.copy;
                    tmpPoint.v1=aManifoldSlice(jj).point(ii).v1.copy;
                    tmpPoint.v2=aManifoldSlice(jj).point(ii).v2.copy;
                    tmpPoint.g1=aManifoldSlice(jj).point(ii).g1.copy;
                    tmpPoint.g2=aManifoldSlice(jj).point(ii).g2.copy;
                    tmpPoint.u=aManifoldSlice(jj).point(ii).u.copy;
                    tmpPoint.r=aManifoldSlice(jj).point(ii).r.copy;
                    tmpPoint.l=aManifoldSlice(jj).point(ii).l.copy;
                    
                    tmpPoint.exitflagNV = -Inf;
                    tmpPoint.exitflagManifold = -Inf;
                    
                    %% disturb manifold point
                    tmpPoint.alpha.values = tmpPoint.alpha.values + aManifoldSlice(jj).stepLength*tngntSpc;
                    
                    
                    
                    %% linear extrapolation of other variables
                    if length(aManifoldSlice(jj).point) > 2
                        tmpPoint.x.values = 2*aManifoldSlice(jj).point(ii).x.values - aManifoldSlice(jj).point(ii - 1).x.values;
                        %                         tmpPoint.alpha.values = 2*aManifoldSlice(jj).point(ii).alpha.values - aManifoldSlice(jj).point(ii - 1).alpha.values;
                        tmpPoint.p.values = 2*aManifoldSlice(jj).point(ii).p.values - aManifoldSlice(jj).point(ii - 1).p.values;
                        tmpPoint.omega.values = 2*aManifoldSlice(jj).point(ii).omega.values - aManifoldSlice(jj).point(ii - 1).omega.values;
                        tmpPoint.w1.values = 2*aManifoldSlice(jj).point(ii).w1.values - aManifoldSlice(jj).point(ii - 1).w1.values ;
                        tmpPoint.w2.values = 2*aManifoldSlice(jj).point(ii).w2.values - aManifoldSlice(jj).point(ii - 1).w2.values;
                    end
                    
                    aManifoldSlice(jj).point(ii+1) = tmpPoint;
                    
                    %% correction
                    
                    offset = aManifoldSlice(jj).point(ii+1).x.index(1)-1;
                    otherVariables=zeros(offset,1);
                    otherVariables(aManifoldSlice(jj).point(ii+1).p.index) = aManifoldSlice(jj).point(ii+1).p.values;
                    
                    
                    rhsAug = @(y)aManifoldSlice(jj).eqAugSys.conFun([otherVariables;...
                        y(aManifoldSlice(jj).point(ii+1).x.index - offset );...
                        aManifoldSlice(jj).point(ii+1).alpha.values + y(aManifoldSlice(jj).point(ii+1).x.nVar+1)*nrmlSpc;... % search direction
                        y(aManifoldSlice(jj).point(ii+1).omega.index - offset - (2+1) + 0*aManifoldSlice(jj).point(ii+1).alpha.nVar);...
                        y(aManifoldSlice(jj).point(ii+1).w1.index - offset - aManifoldSlice(jj).point(ii+1).alpha.nVar + 2 - 1);% + 0*(2+1) + 0*aManifoldSlice(jj).point(ii+1).alpha.nVar);...
                        y(aManifoldSlice(jj).point(ii+1).w2.index - offset - aManifoldSlice(jj).point(ii+1).alpha.nVar + 2 - 1)]);% + 0*(2+1) + 0*aManifoldSlice(jj).point(ii+1).alpha.nVar)]);
                    
                    
                    
                    x0 = [aManifoldSlice(jj).point(ii+1).x.values;
                        0;
                        aManifoldSlice(jj).point(ii+1).omega.values;
                        aManifoldSlice(jj).point(ii+1).w1.values;
                        aManifoldSlice(jj).point(ii+1).w2.values];
                    
                    
                    options = optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',200000,'display','off','TolFun',1e-9,'TolX',1e-6);
                    
                    
                    [x,~,exitflag,output] = fsolve(rhsAug,x0,options);
                    
                    if (output.iterations < 200) && (exitflag>0) && (x(aManifoldSlice(jj).point(ii+1).x.nVar+1)<aManifoldSlice(jj).maxStepLength)
                        aManifoldSlice(jj).stepLength=min(aManifoldSlice(jj).stepLength*1.2,aManifoldSlice(jj).maxStepLength);
                    else
                        aManifoldSlice(jj).stepLength=aManifoldSlice(jj).stepLength/2;
                    end
                    
                    if aManifoldSlice.debugFlag && (exitflag < 1)
                        warning('exitflag was %d during correction', exitflag)
                    end
                    aManifoldSlice(jj).point(ii+1).x.values = x(aManifoldSlice(jj).point(ii+1).x.index - offset );
                    aManifoldSlice(jj).point(ii+1).alpha.values = aManifoldSlice(jj).point(ii+1).alpha.values + x(aManifoldSlice(jj).point(ii+1).x.nVar+1)*aManifoldSlice(jj).point(ii+1).r.values;
                    aManifoldSlice(jj).point(ii+1).omega.values = x(aManifoldSlice(jj).point(ii+1).omega.index - offset - (2+1) + aManifoldSlice(jj).point(ii+1).alpha.nVar);
                    aManifoldSlice(jj).point(ii+1).w1.values = x(aManifoldSlice(jj).point(ii+1).w1.index - offset - aManifoldSlice(jj).point(ii+1).alpha.nVar + 2 - 1);
                    aManifoldSlice(jj).point(ii+1).w2.values = x(aManifoldSlice(jj).point(ii+1).w2.index - offset - aManifoldSlice(jj).point(ii+1).alpha.nVar + 2 - 1);
                    aManifoldSlice(jj).point(ii+1).exitflagManifold = exitflag;
                    
                    %% find normal vector
                    
                    offset=aManifoldSlice(jj).point(ii+1).x.index(1)-1;
                    
                    otherVariables=zeros(offset,1);
                    otherVariables(aManifoldSlice(jj).point(ii+1).p.index)=aManifoldSlice(jj).point(ii+1).p.values;
                    
                    maniPoint = [aManifoldSlice(jj).point(ii+1).x.values;...
                        aManifoldSlice(jj).point(ii+1).alpha.values;...
                        aManifoldSlice(jj).point(ii+1).omega.values;...
                        aManifoldSlice(jj).point(ii+1).w1.values;...
                        aManifoldSlice(jj).point(ii+1).w2.values];
                    
                    
                    
                    x0 = [aManifoldSlice(jj).point(ii+1).v1.values;
                        aManifoldSlice(jj).point(ii+1).v2.values;
                        aManifoldSlice(jj).point(ii+1).g1.values;
                        aManifoldSlice(jj).point(ii+1).g2.values;
                        aManifoldSlice(jj).point(ii+1).u.values;
                        aManifoldSlice(jj).point(ii+1).r.values];
                    
                    rhsNV = @(y)aManifoldSlice(jj).eqNVSys.conFun([otherVariables; maniPoint; y]);
                    
                    
                    options = optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxIter',1000,'MaxFunEvals',2000000,'display','off','TolFun',1e-9,'TolX',1e-6);
                    
                    [x,~,exitflag] = fsolve(rhsNV,x0,options);
                    
                    if aManifoldSlice.debugFlag && (exitflag < 1)
                        warning('exitflag was %d during NV search', exitflag)
                    end
                    
                    aManifoldSlice(jj).point(ii+1).v1.values = x(aManifoldSlice(jj).point(ii+1).v1.index - offset - length(maniPoint));
                    aManifoldSlice(jj).point(ii+1).v2.values = x(aManifoldSlice(jj).point(ii+1).v2.index - offset - length(maniPoint));
                    aManifoldSlice(jj).point(ii+1).g1.values = x(aManifoldSlice(jj).point(ii+1).g1.index - offset - length(maniPoint));
                    aManifoldSlice(jj).point(ii+1).g2.values = x(aManifoldSlice(jj).point(ii+1).g2.index - offset - length(maniPoint));
                    aManifoldSlice(jj).point(ii+1).u.values = x(aManifoldSlice(jj).point(ii+1).u.index - offset - length(maniPoint));
                    aManifoldSlice(jj).point(ii+1).r.values = x(aManifoldSlice(jj).point(ii+1).r.index- offset - length(maniPoint));
                    aManifoldSlice(jj).point(ii+1).exitflagNV = exitflag;
                    
                    if aManifoldSlice(jj).showStepsFlag
                        plot(aManifoldSlice(jj).point(ii+1).alpha.values(aManifoldSlice(jj).freeParamIndices(1)),aManifoldSlice(jj).point(ii+1).alpha.values(aManifoldSlice(jj).freeParamIndices(2)),'bo')
                        quiver(aManifoldSlice(jj).point(ii+1).alpha.values(aManifoldSlice(jj).freeParamIndices(1)),aManifoldSlice(jj).point(ii+1).alpha.values(aManifoldSlice(jj).freeParamIndices(2)),aManifoldSlice(jj).point(ii+1).r.values(aManifoldSlice(jj).freeParamIndices(1)),aManifoldSlice(jj).point(ii+1).r.values(aManifoldSlice(jj).freeParamIndices(2)),'r')
                        quiver(aManifoldSlice(jj).point(ii+1).alpha.values(aManifoldSlice(jj).freeParamIndices(1)),aManifoldSlice(jj).point(ii+1).alpha.values(aManifoldSlice(jj).freeParamIndices(2)),tngntSpc(aManifoldSlice(jj).freeParamIndices(1)),tngntSpc(aManifoldSlice(jj).freeParamIndices(2)),'m');
                        drawnow
                    end
                    
                    if ~all(aManifoldSlice(jj).point(ii+1).alpha.values <= aManifoldSlice(jj).upperBoxCons)...
                            || ~all(aManifoldSlice(jj).point(ii+1).alpha.values >= aManifoldSlice(jj).lowerBoxCons)
                        disp('boundary hit')
                        disp([aManifoldSlice(jj).lowerBoxCons,aManifoldSlice(jj).point(ii+1).alpha.values,aManifoldSlice(jj).upperBoxCons])
                        
                        break
                    end
                end
            end
        end
        
        
    % ======================================================================
    %> @brief plots results of numerical continuation of critical manifolds
    %>
    %> @param aManifoldSlice instance of this class
    %>
    %> @return handle vector of plot handles allowing later manipulation of plots
    % ======================================================================
        
        function handle = plot(aManifoldSlice)
            for jj = 1:length(aManifoldSlice)
                hold on
                
                alpha=NaN(2,length(aManifoldSlice(jj).point));
                
                for ii = 1:length(aManifoldSlice(jj).point)
                    if aManifoldSlice(jj).point(ii).exitflagManifold>0
                        alpha(1,ii) = aManifoldSlice(jj).point(ii).alpha.values(aManifoldSlice(jj).freeParamIndices(1));
                        alpha(2,ii) = aManifoldSlice(jj).point(ii).alpha.values(aManifoldSlice(jj).freeParamIndices(2));
                    end
                end
                
                xlabel(aManifoldSlice(jj).point(1).alpha.names{aManifoldSlice(jj).freeParamIndices(1)})
                ylabel(aManifoldSlice(jj).point(1).alpha.names{aManifoldSlice(jj).freeParamIndices(2)})
                

                handle(jj) = plot(alpha(1,:), alpha(2,:), '-');
            end
        end
        
        function quiver(aManifoldSlice)
            for jj = 1:length(aManifoldSlice)
                hold on
                
                alpha=NaN(2,length(aManifoldSlice(jj).point));
                r=NaN(2,length(aManifoldSlice(jj).point));
                for ii = 1:length(aManifoldSlice(jj).point)
                    if aManifoldSlice(jj).point(ii).exitflagNV>0
                        alpha(1,ii) = aManifoldSlice(jj).point(ii).alpha.values(aManifoldSlice(jj).freeParamIndices(1));
                        alpha(2,ii) = aManifoldSlice(jj).point(ii).alpha.values(aManifoldSlice(jj).freeParamIndices(2));
                        r(1,ii) = aManifoldSlice(jj).point(ii).r.values(aManifoldSlice(jj).freeParamIndices(1));
                        r(2,ii) = aManifoldSlice(jj).point(ii).r.values(aManifoldSlice(jj).freeParamIndices(2));
                        r(:,ii) = 0.1*r(:,ii)/norm(r(:,ii),2);
                    end
                end
                quiver(alpha(1,:), alpha(2,:),r(1,:),r(2,:),'b-')
            end
        end
    end
end

