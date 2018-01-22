%> @file  varCollection.m
%> @brief creates variable collections for different manifold types
%> @author Jonas Otten
%> @date 18 Jul 2017
%======================================================================

%======================================================================
%> @brief creates variable collections for different manifold types
%>
%> the structs are allways created the same, but for some manifolds 
%> with empty entries
%>
%> @param type requested manifold type
%> @param offset index offset
%> @param x state vector, instance of class VariableVector
%> @param alpha uncertain variables, instance of class VariableVector
%> @param p certain optimization variables, instance of class VariableVector
%>
%> @retval vars variables vector to return
%======================================================================

function [ vars ] = varCollection( type, offset, x, alpha, p)
%VARCOLLETION creates variable collections for different manifold types


vars.x=VariableVector(x.values,offset,x.names);
vars.alpha=VariableVector(alpha.values,offset+x.nVar,alpha.names);
vars.p=VariableVector(p.values,offset+x.nVar+alpha.nVar,p.names);
vars.omega=VariableVector([], offset+x.nVar,{});
vars.w1=VariableVector([], offset+x.nVar,{});
vars.w2=VariableVector([], offset+x.nVar,{});
vars.v1=VariableVector([], offset+x.nVar,{});
vars.v2=VariableVector([], offset+x.nVar,{});
vars.g1=VariableVector([], offset+x.nVar,{});
vars.g2=VariableVector([], offset+x.nVar,{});
vars.u=VariableVector([], offset+x.nVar,{});
vars.r=VariableVector([],offset+x.nVar,{});
vars.l=VariableVector([],offset+x.nVar,{});
vars.maxEig=[];
vars.eigs=[];

switch type
    case 'nominal'
        vars.nVar=x.nVar+alpha.nVar;    
    case 'fold'
        vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar+p.nVar,{'w1'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+2*x.nVar+alpha.nVar+p.nVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+2*x.nVar+2*alpha.nVar+p.nVar,{'l'});
        vars.nVar=2*x.nVar+2*alpha.nVar+p.nVar+1;
%     case 'modfold'
%         vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar,{'w1'});
%         vars.v1=VariableVector(zeros(x.nVar,1),offset+2*x.nVar+alpha.nVar,{'v1'});
%         vars.g1=VariableVector(zeros(1,1),offset+3*x.nVar+alpha.nVar,{'g1'});
%         vars.u=VariableVector(zeros(x.nVar,1),offset+3*x.nVar+alpha.nVar+1,{'u'});
%         vars.r=VariableVector(zeros(alpha.nVar,1),offset+4*x.nVar+alpha.nVar+1,{'r'});
%         vars.l=VariableVector(zeros(1,1),offset+4*x.nVar+2*alpha.nVar+1,{'l'});
%         vars.nVar=4*x.nVar+2*alpha.nVar+2;

%HIER MODIFIZIERT UM AUCH DIE INDIZES DER FOLGENDEN VARIABLEN ANZUPASSEN
    case 'modfold'
        vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar+p.nVar,{'w1'});
        vars.v1=VariableVector(zeros(x.nVar,1),offset+2*x.nVar+alpha.nVar+p.nVar,{'v1'});
        vars.g1=VariableVector(zeros(1,1),offset+3*x.nVar+alpha.nVar+p.nVar,{'g1'});
        vars.u=VariableVector(zeros(x.nVar,1),offset+3*x.nVar+alpha.nVar+1+p.nVar,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+4*x.nVar+alpha.nVar+1+p.nVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+4*x.nVar+2*alpha.nVar+1+p.nVar,{'l'});
        vars.nVar=4*x.nVar+2*alpha.nVar+2+p.nVar;
%Bis hier!
        
    case 'hopf'
        vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar+p.nVar,{'omega'});
        vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+p.nVar+1,{'w1'});
        vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+p.nVar+1,{'w2'});
        vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+p.nVar+1,{'v1'});
        vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+p.nVar+1,{'v2'});
        vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+p.nVar+1,{'g1'});
        vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+p.nVar+2,{'g2'});
        vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+p.nVar+3,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+p.nVar+3,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+p.nVar+3,{'l'});
        vars.nVar=6*x.nVar+2*alpha.nVar+p.nVar+4;
%     case 'modhopf'
%         vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar,{'omega'});
%         vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+1,{'w1'});
%         vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+1,{'w2'});
%         vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+1,{'v1'});
%         vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+1,{'v2'});
%         vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+1,{'g1'});
%         vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+2,{'g2'});
%         vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+3,{'u'});
%         vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+3,{'r'});
%         vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+3,{'l'});
%         vars.nVar=6*x.nVar+2*alpha.nVar+4;
        
% modifizieren, um auch hier das p mit zu berücksichtigen

        case 'modhopf'
        vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar+p.nVar,{'omega'});
        vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+1+p.nVar,{'w1'});
        vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+1+p.nVar,{'w2'});
        vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+1+p.nVar,{'v1'});
        vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+1+p.nVar,{'v2'});
        vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+1+p.nVar,{'g1'});
        vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+2+p.nVar,{'g2'});
        vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+3+p.nVar,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+3+p.nVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+3+p.nVar,{'l'});
        vars.nVar=6*x.nVar+2*alpha.nVar+4+p.nVar;
    otherwise
        error('unknown type to generate variable list')
end


end

