classdef EqualityConstraint < handle
    %EQUALITYCONSTRAINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
     end
    
    properties(SetAccess=protected)
        conFun % function handle
        eqIndex
        nEqs
        status=0;
    end
    
    methods
        function anEqCon=EqualityConstraint(conFunHandle,nEqs,vars,eqOffset)
            % constructor
            % fill the properties
            anEqCon.conFun=conFunHandle;
            anEqCon.vars=vars;
            anEqCon.eqIndex=eqOffset+(1:nEqs);
            anEqCon.nEqs=nEqs;
        end
        
        function anEqCon=shiftIndex(anEqCon,eqShift,varShift)
            % shift indices if other constraints have been removed            

            % do this at first for the equation indices
            lowestEqIndex = min(anEqCon.eqIndex);
            if ((eqShift + lowestEqIndex) < 1)
                error('shifting the equation index would lead to negative indices');
            end
            anEqCon.eqIndex=anEqCon.eqIndex + eqShift;
            
            % and use the lower rank methode for the variables
            anEqCon.vars.shiftIndex(varShift);
            
        end        
    end    
end

