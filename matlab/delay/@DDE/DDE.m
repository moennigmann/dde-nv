classdef DDE < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % handles to manifold defining functions
        foldManiHandle = []
        modfoldManiHandle = []
        hopfManiHandle = []
        modhopfManiHandle = []

        % handles to NV defining functions
        foldNVHandle = []
        modfoldNVHandle = []
        hopfNVHandle = []
        modhopfNVHandle = []
            
    end

    properties(SetAccess=protected)
                        
        rhs % handle for right-hand-side of DDE
        delays % handle of function for delays
        
        uncParam % uncertain parameters: object of class VariableVector
        certOptParam % certain optimization parameters: object of class VariableVector
                
        ntau % number of delays
        
    end
    
    methods     
        function aDDE = DDE(rhsHandle,delays,xNomGuess,uncParam,certOptParam) % constructor
            
            if nargin < 5
                certOptParam = VariableVector([],xNomGuess.nVar+uncParam.nVar,[]);
            end
            
            %% write inputs into object properties
            
            % starting with constant properties
            aDDE.rhs = rhsHandle;
            
            aDDE.delays = delays; 
            
            aDDE.ntau = length(delays(xNomGuess.values,uncParam.values));
              
            
            %% define state specific properties
            
            
            %% define parameter specific properties
            aDDE.uncParam = uncParam;
            aDDE.certOptParam = certOptParam;
            
        end
        
        function [handleMani, handleNV] = getHandles(aDDE,type)
            switch type
                case 'fold'
                    handleMani = aDDE.foldManiHandle;
                    handleNV  = aDDE.foldNVHandle;
                case 'modfold'
                    handleMani = aDDE.modfoldManiHandle;
                    handleNV = aDDE.modfoldNVHandle;
                case 'hopf'
                    handleMani = aDDE.hopfManiHandle;
                    handleNV = aDDE.hopfNVHandle;
                case 'modhopf'
                    handleMani = aDDE.modhopfManiHandle;
                    handleNV = aDDE.modhopfNVHandle;
                otherwise
                    error('unknown manifold type requested');
            end
            if all(size(handleMani) == [0 0])
                error('the requested manifold function handle (%s) was not found',type)
            end
            if all(size(handleMani) == [0 0])
                error('the requested normal vector system handle was not found')
            end
        end
        
    end
    
end

