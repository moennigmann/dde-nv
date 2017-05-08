function [ vars ] = varCollection( type, offset, x, alpha, p)
%VARCOLLETION creates variable collections for different manifold types


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
vars.u=VariableVector([], offset+x.nVar,{});
vars.r=VariableVector([],offset+x.nVar,{});
vars.l=VariableVector([],offset+x.nVar,{});

switch type
    case 'nominal'
        vars.nVar=x.nVar+alpha.nVar;    
    case 'fold'
        vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar,{'w1'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+2*x.nVar+alpha.nVar,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+2*x.nVar+2*alpha.nVar,{'l'});
        vars.nVar=2*x.nVar+2*alpha.nVar+1;
    case 'modfold'
        vars.w1=VariableVector(ones(x.nVar,1),offset+x.nVar+alpha.nVar,{'w1'});
        vars.v1=VariableVector(zeros(x.nVar,1),offset+2*x.nVar+alpha.nVar,{'v1'});
        vars.g1=VariableVector(zeros(1,1),offset+3*x.nVar+alpha.nVar,{'g1'});
        vars.u=VariableVector(zeros(x.nVar,1),offset+3*x.nVar+alpha.nVar+1,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+4*x.nVar+alpha.nVar+1,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+4*x.nVar+2*alpha.nVar+1,{'l'});
        vars.nVar=4*x.nVar+2*alpha.nVar+2;
    case 'hopf'
        vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar,{'omega'});
        vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+1,{'w1'});
        vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+1,{'w2'});
        vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+1,{'v1'});
        vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+1,{'v2'});
        vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+1,{'g1'});
        vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+2,{'g2'});
        vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+3,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+3,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+3,{'l'});
        vars.nVar=6*x.nVar+2*alpha.nVar+4;
    case 'modhopf'
        vars.omega=VariableVector(zeros(1,1), offset+x.nVar+alpha.nVar,{'omega'});
        vars.w1=VariableVector(ones(x.nVar,1), offset+x.nVar+alpha.nVar+1,{'w1'});
        vars.w2=VariableVector(zeros(x.nVar,1), offset+2*x.nVar+alpha.nVar+1,{'w2'});
        vars.v1=VariableVector(zeros(x.nVar,1), offset+3*x.nVar+alpha.nVar+1,{'v1'});
        vars.v2=VariableVector(zeros(x.nVar,1), offset+4*x.nVar+alpha.nVar+1,{'v2'});
        vars.g1=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+1,{'g1'});
        vars.g2=VariableVector(zeros(1,1), offset+5*x.nVar+alpha.nVar+2,{'g2'});
        vars.u=VariableVector(zeros(x.nVar,1), offset+5*x.nVar+alpha.nVar+3,{'u'});
        vars.r=VariableVector(zeros(alpha.nVar,1),offset+6*x.nVar+alpha.nVar+3,{'r'});
        vars.l=VariableVector(zeros(1,1),offset+6*x.nVar+2*alpha.nVar+3,{'l'});
        vars.nVar=6*x.nVar+2*alpha.nVar+4;
    otherwise
        error('unknown type to generate variable list')
end


end

