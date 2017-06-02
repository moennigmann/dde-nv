classdef VariableVector < handle
    % The instances of this class contain information for a variable used
    % in DDENLP
    %   Those instances collect human-readable information on the variables
    %   (their names and valueus) as well as the information which is used
    %   by other objects to adress those variables.

    
    
    properties
              values % column vector containing the numerical values of the Variables
    end
    
    properties (SetAccess = protected)
        names % column array of strings containing the description of the 
                % variables
        index % indexes of variable within superordinate concatenation of
                % variables
        nVar% number of variables
    end
    
    methods
        function aVariables = VariableVector( values,offset,nameInput ) % constructor
            % This function constructs instances of the class Variable Vector
            %
            
            % check if 'values' is a column vector. If yes, assign 'nVar'
            if size( values,2 ) <= 1
                aVariables.nVar = size( values,1 );
            else
                error('''values'' (numerical values of variables) has to be a column vector')
            end
            
            % check if nameInput was an argument of the constructor. If
            % not, create empty name array
            if nargin < 3
                nameInput=repmat({'not assigned'},aVariables.nVar,1);
            else
                if all( size( nameInput ) == [1,1]) && ~all( size( values ) == [1,1] )
                    name=nameInput{1};
                    for ii = 1:aVariables.nVar
                        nameInput{ii,1} = {[name,'_',num2str(ii)]};
                    end
                end
            end
            
            % check if nameInput is a column vectors
            if size(nameInput,2) > 1
                    error('''nameInput'' has to be a column vector')
            end    
            
            % check if number of names and values are identical
            if size(nameInput,1) == size(values,1)
                aVariables.names = nameInput;
                aVariables.index = offset+1:offset+aVariables.nVar;
                aVariables.values = values;
            else
                error('length of ''nameInput'' and ''values'' have to be identical')
            end
        end
        
        function shiftIndex(aVariableVec , shift)
            % this method shifts the index of the Variables in the
            % instance of VariableVector 
            lowestVarIndex = min(aVariableVec.index);
            if ((shift + lowestVarIndex) < 1)
                error('shifting the variable index would lead to negative indices');
            end
            aVariableVec.index = aVariableVec.index + shift;
        end
            
    end
end
    