%> @file EqualityConstraint.m
%> @brief (abstract) class that inherits its properities to equality constraints
% ======================================================================
%> @brief (abstract) class that inherits its properities to equality constraints
%
%>   The properties are necessary to manage mutliple classes of equality constraints
% ======================================================================

classdef EqualityConstraint < handle
    
    properties
        %> collection of variables of class VariableVector on which the instance of EqualityConstraint depends on
        vars
     end
    
    properties(SetAccess=protected)
        %> function handle representin the constraint function
        conFun
        %> indices of the equation in the superordinate constraint optimization problem
        eqIndex
        %> number of equation in this instance of EqualityConstraint
        nEqs
        %> initialization status of this instance of EqualityConstraint
        status=0;
    end
    
    methods
        
                        
    % ======================================================================
    %> @brief Class constructor
    %>
    %> This function constructs instances of the class EqualityConstraint
    %>
    %> @param conFunHandle function handle for equality constraint
    %> @param nEqs number of equations in this equality constraint
    %> @param vars collection of instances of VariableVector
    %> @param eqOffset offset for equation indices
    % 
    %> @return instance of the EqualityConstraint class.
    % ======================================================================
        
        function anEqCon=EqualityConstraint(conFunHandle,nEqs,vars,eqOffset)
            % constructor
            % fill the properties
            anEqCon.conFun=conFunHandle;
            anEqCon.vars=vars;
            anEqCon.eqIndex=eqOffset+(1:nEqs);
            anEqCon.nEqs=nEqs;
        end
        
        
	% ======================================================================
    %> @brief shift the index of the equations stored in this equality constraint
    %>
    %> @param anEqCon Instance of EqualityConstraint, where the variable equation will be shifted
    %> @param eqShift how far the index of equations will be shifted
    %> @param varShift how far the index of variables will be shifted
    %>
    %> @return instance of the EqualityConstraint class.
    % ====================================================================== 
        
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

