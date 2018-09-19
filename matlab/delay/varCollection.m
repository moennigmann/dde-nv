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

function [ vars ] = varCollection( type, offset, x, alpha, p, algVars)
%VARCOLLETION creates variable collections for different manifold types


if algVars
    p = VariableVector(p.values,offset+x.nVar+alpha.nVar,p.names);
    pnVar=p.nVar;
else
    pnVar = 0;
end
    
vars.x=VariableVector(x.values,offset,x.names);
vars.alpha=VariableVector(alpha.values,offset+x.nVar,alpha.names);
vars.p=p;
vars.omega=VariableVector([], offset+x.nVar,{});
vars.w1=VariableVector([], offset+x.nVar,{});
vars.w2=VariableVector([], offset+x.nVar,{});
vars.v1=VariableVector([], offset+x.nVar,{});
vars.v2=VariableVector([], offset+x.nVar,{});
vars.g1=VariableVector([], offset+x.nVar,{});
vars.g2=VariableVector([], offset+x.nVar,{});
vars.k=VariableVector([], offset+x.nVar,{});
vars.u=VariableVector([], offset+x.nVar,{});
vars.r=VariableVector([],offset+x.nVar,{});
vars.l=VariableVector([],offset+x.nVar,{});
vars.maxEig=[];
vars.eigs=[];

switch type
    case 'nominal'
        vars.nVar=x.nVar+alpha.nVar;    
    case 'fold'
        vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar+pnVar,{'w1'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+2*x.nVar+alpha.nVar+pnVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+2*x.nVar+2*alpha.nVar+pnVar,{'l'});
        vars.nVar=2*x.nVar+2*alpha.nVar+pnVar+1;
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
        vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar+pnVar,{'w1'});
        vars.v1=VariableVector(zeros(x.nVar,1),offset+2*x.nVar+alpha.nVar+pnVar,{'v1'});
        vars.g1=VariableVector(zeros(1,1),offset+3*x.nVar+alpha.nVar+pnVar,{'g1'});
        vars.k=VariableVector(zeros(pnVar,1),offset+3*x.nVar+alpha.nVar+1+pnVar,{});
        vars.u=VariableVector(zeros(x.nVar,1),offset+3*x.nVar+alpha.nVar+1+2*pnVar,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+4*x.nVar+alpha.nVar+1+2*pnVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+4*x.nVar+2*alpha.nVar+1+2*pnVar,{'l'});
        vars.nVar=4*x.nVar+2*alpha.nVar+2+2*pnVar;

        
    case 'hopf'
        vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar+pnVar,{'omega'});
        vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+pnVar+1,{'w1'});
        vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+pnVar+1,{'w2'});
        vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+pnVar+1,{'v1'});
        vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+pnVar+1,{'v2'});
        vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+pnVar+1,{'g1'});
        vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+pnVar+2,{'g2'});
        vars.k=VariableVector(zeros(pnVar,1), offset+5*x.nVar+alpha.nVar+pnVar+3,{'k'});
        vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+2*pnVar+3,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+2*pnVar+3,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+2*pnVar+3,{'l'});
        vars.nVar=6*x.nVar+2*alpha.nVar+2*pnVar+4;
%Bis hier!
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
        vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar+pnVar,{'omega'});
        vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+1+pnVar,{'w1'});
        vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+1+pnVar,{'w2'});
        vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+1+pnVar,{'v1'});
        vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+1+pnVar,{'v2'});
        vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+1+pnVar,{'g1'});
        vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+2+pnVar,{'g2'});
        vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+3+pnVar,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+3+pnVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+3+pnVar,{'l'});
        vars.nVar=6*x.nVar+2*alpha.nVar+4+pnVar;
    otherwise
        error('unknown type to generate variable list')
end


end

