classdef StStConstraint < EqualityConstraint
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(SetAccess=protected)

    end
    
    methods
        function aStStCon=StStConstraint(aDDE,vars)
            
            conEq=@(y)aDDE.rhs(y(vars.x.index),...
                repmat(y(vars.x.index),1,aDDE.ntau),...
                y(vars.alpha.index),...
                y(vars.p.index));
            aStStCon=aStStCon@EqualityConstraint(conEq,vars.x.nVar,vars,0);
        end
        
        
        function aStStCon=initStStConstraint(aStStCon,options)
            %create function handle to search for initial state while
            %holding the parameters alpha fixed
            
            indexX=aStStCon.vars.x.index;
            alpha=aStStCon.vars.alpha.values;
            p=aStStCon.vars.p.values;
            conEqFixedAlpha=@(y)aStStCon.conFun([y(indexX);alpha;p]);
            
            x0=aStStCon.vars.x.values;
            % try to initialize steady state constraint
            [x,~,exitflag]=fsolve(conEqFixedAlpha,x0,options);
            aStStCon.vars.x.values=x;
            if exitflag>0
                aStStCon.status=1;
                fprintf('steady state constraint for nominal point successfully initialized!\n')
            else
                warning(['initialization of steady state contsraint with default options for fsolve not successful, exitflag = ',num2str(exitflag)]);
            end
        end
    end
    
end

