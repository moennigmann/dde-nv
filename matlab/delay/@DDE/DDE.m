%> @file DDE.m
%> @brief The instances of this class contain information for a DDE used in DDENLP
% ======================================================================
%> @brief The instances of this class contain information for a DDE used in DDENLP
%
%>   Those instances collect the relevant data of the DDE describing the system to be optimized
% ======================================================================


classdef DDE < handle
    
    properties
        % handles to manifold defining functions
        %> function handle for fold bifurcation manifold
        foldManiHandle = []
        %> function handle for modified fold bifurcation manifold
        modfoldManiHandle = []        
        %> function handle for hopf bifurcation manifold
        hopfManiHandle = []
        %> function handle for modified hopf bifurcation manifold
        modhopfManiHandle = []

        % handles to NV defining functions
        %> function handle for fold bifurcation normal vector system
        foldNVHandle = []
        %> function handle for modified fold bifurcation normal vector system
        modfoldNVHandle = []        
        %> function handle for hopf bifurcation normal vector system
        hopfNVHandle = []
        %> function handle for modified hopf bifurcation normal vector system
        modhopfNVHandle = []
            
    end

    properties(SetAccess=protected)
        
        %> function handle for right-hand-side of DDE              
        rhs
        %> function handle of function for delays
        delays
        
        %> uncertain parameters: object of class VariableVector
        uncParam 
        %> certain optimization parameters: object of class VariableVector
        certOptParam 
        %> number of delays
        ntau 
        
    end
    
    methods
                
    % ======================================================================
    %> @brief Class constructor
    %>
    %> This function constructs instances of the class DDE
    %> @param rhsHandle function handle for right hand side of DDE
    %> @param delays function handle for delays
    %> @param xNomGuess Guess for a nominal steady state vector of class VariableVector
    %> @param uncParam uncertain nominal paraemeter vector of class VariableVector
    %> @param certOptParam parameter vector with certain optimization variables of class VariableVector
    %>
    %> @return instance of the DDE class.
    % ======================================================================
        
        
        function aDDE = DDE(rhsHandle,delays,xNomGuess,uncParam,certOptParam) % constructor
            
            if nargin < 5
                certOptParam = VariableVector([],xNomGuess.nVar+uncParam.nVar,[]);
            end
            
            %% write inputs into object properties
            
            % starting with constant properties
            aDDE.rhs = rhsHandle;
            
            aDDE.delays = delays; 
            
            aDDE.ntau = length(delays(xNomGuess.values,uncParam.values,certOptParam.values));
              
            
            %% define state specific properties
            
            
            %% define parameter specific properties
            aDDE.uncParam = uncParam;
            aDDE.certOptParam = certOptParam;
            
        end
        
    % ======================================================================
    %> @brief hands back function handles of a DDE instance
    %>
    %> Those function handles come in handy for the formulation of normal vector constraints 
    %> @param aDDE function handle for right hand side of DDE
    %> @param type function handle for delays
    %>
    %> @return handleMani function handle with function g(.)=0 describing a critical manifold
    %> @return handleNV function handle with function h(.)=0 describing normal vectors on critical manifold
    % ======================================================================
            
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
            if all(size(handleNV) == [0 0])
                error('the requested normal vector system handle was not found')
            end
        end
        
    end
    
end

