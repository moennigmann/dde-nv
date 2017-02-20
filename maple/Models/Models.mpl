Models := module () export Aiello1992, Agrawal1982, Disney2002, Disney2003, 
Goswami1996, John1994, Lehman1994, Luzyanina2001, Minimal, Olsen1983, 
Scott1990, SOFC, Teymour1989b; Aiello1992 := module () export Sys; Sys := 
module () local aSys; export getSys; aSys := table([("ExplicitAEs")=[gamma = 1
],("Parameters")=[alpha = 1, beta = 1],("Delays")=[tau[1] = 2-exp(-x2)],(
"AlgVars")=[],("DynVars")=[x1, x2],("ODEs")=[`x1'` = alpha*x2-gamma*x1-alpha*
exp(gamma*tau[1])*x2Del, `x2'` = alpha*exp(gamma*tau[1])*x2Del-beta*x2^2],(
"DelVars")=[x1Del, x2Del],("AEs")=[]]); getSys := proc () return eval(aSys) 
end proc; end module; end module; Agrawal1982 := module () export CSTR, 
DimensionalCSTR; CSTR := module () local Sys; export getSys, getSolutions; Sys
:= table([("ExplicitAEs")=[k = 1, K = .12, a = 5.4, b = 180, M = -(x2-1)*exp(
SF*x2/K), Sigma = M*(a+b*SF)/(a-b*SF*x2+b*SF)],("Parameters")=[Da = .910365, 
SF = .3],("AlgVars")=[],("DynVars")=[x1, x2],("ODEs")=[`x1'` = -x1+M*Da*x1, 
`x2'` = -x2+Sigma*Da*x1],("AEs")=[]]); getSys := proc () return eval(Sys) end
proc; getSolutions := proc () local Sol1, solutions; Sol1 := [x1 = .17747, x2
= .87753, Da = .910365, SF = .3, k = 1., K = .12, a = 5.4, b = 180., M = 1.098\
489246, Sigma = 5.431465682]; solutions := [Sol1]; return solutions end proc;
end module; DimensionalCSTR := module () local Sys; export getSys, 
getSolutions; Sys := table([("ExplicitAEs")=[k = 1, K = .12, a = 5.4, b = 180,
c1 = 5, V = 51.3*Pi, mu = k*S*exp(-S/K), sigma = mu/(a+b*S), c2 = -70.3+270.3*
SF],("Parameters")=[F = 3.0, SF = .3],("AlgVars")=[],("DynVars")=[X, S],(
"ODEs")=[`X'` = -F/V*X+mu*X, `S'` = F/V*(SF-S)-sigma*X],("AEs")=[]]); getSys 
:= proc () return eval(Sys) end proc; getSolutions := proc () local Sol1, Sol2
, solutions; Sol1 := [F = 3, SF = .3, X = 2.6201, S = .22443e-1, k = 1., K = .\
12, a = 5.4, b = 180., c1 = 5., V = 161.1637032, mu = .1861474040e-1, sigma =
.1971954779e-2, c2 = 10.79]; Sol2 := [X = 4.22744, S = .7385282e-1, F = 6.4320\
89, SF = .3, k = 1., K = .12, a = 5.4, b = 180., c1 = 5., V = 161.1637032, mu
= .3991028368e-1, sigma = .2134981007e-2, c2 = 10.79, w1_Hopf[1] = -.\
1102549550e-7, w1_Hopf[2] = -.4923686503e-1, w2_Hopf[1] = .9987871301, w2_Hopf
[2] = -.2294538897e-6, omega_Hopf = .4330884089e-1]; solutions := [Sol1, Sol2]
; return solutions end proc; end module; end module; Disney2002 := module () 
export Sys; Sys := module () local Syst; export getSys; Syst := table([(
"ExplicitAlgEqns")=[DINV = 0, WIP = ORATE+ORATEaux1+ORATEaux2+ORATEaux3, 
VCON_t = CONSpar+VCON+1/(1+Tq)*(G*CONSpar-VCON)-VCONaux-1/(1+Tq)*(G*CONSaux-
VCONaux), AVCON_t = AVCON+1/(1+Ta)*(VCON_t-AVCON), DWIPaux = AVCON*Tpbar, 
COMRATE_t = ORATEaux3, AINV_t = AINV+COMRATE_t-CONSpar, EINVaux = DINV-AINV, 
ORATE_t = AVCON+EINVaux/Ti+(DWIPaux-WIP)/Tw],("Parameters")=[CONSpar = 0, Ta =
8, G = 2, Tq = 4, Tpbar = 4, Ti = 4, Tw = 4],("AlgEqns")=[],("AlgVars")=[],(
"DynVars")=[VCON, AVCON, CONSt, CONSaux, VCONaux, ORATEaux1, ORATEaux2, 
ORATEaux3, COMRATE, ORATE, AINV, R],("DynEqns")=[VCON = VCON_t, AVCON = 
AVCON_t, CONSt = CONSpar, CONSaux = CONSpar, VCONaux = VCON, ORATEaux1 = ORATE
, ORATEaux2 = ORATEaux1, ORATEaux3 = ORATEaux2, COMRATE = COMRATE_t, ORATE = 
ORATE_t, AINV = AINV_t, R = R+1/(1+Tq)*(G*CONSpar-R)]]); getSys := proc () 
return eval(Syst) end proc; end module; end module; Disney2003 := module () 
export Sys; Sys := module () local Syst; export getSys; Syst := table([(
"ExplicitAlgEqns")=[DINV = 0, WIP = ORATE+ORATEaux1, AVCON_t = AVCON+(CONSpar-
AVCON)/(1+Ta), DWIPaux = AVCON*Tpbar, COMRATE_t = ORATEaux1, AINV_t = AINV+
COMRATE_t-CONSpar, EINVaux = DINV-AINV, ORATE_t = AVCON+EINVaux/Ti+(DWIPaux-
WIP)/Tw, DINVD = 0, WIPD = ORATED+ORATEaux1D+ORATEaux2D, AVCON_tD = AVCOND+(
ORATE_t-AVCOND)/(1+Ta2), DWIPauxD = AVCOND*Tpbar2, COMRATE_tD = ORATEaux2D, 
AINV_tD = AINVD+COMRATE_tD-ORATE_t, EINVauxD = DINVD-AINVD, ORATE_tD = AVCOND+
EINVauxD/Ti2+(DWIPauxD-WIPD)/Tw2, DINVM = 0, WIPM = ORATEM+ORATEaux1M+
ORATEaux2M+ORATEaux3M, AVCON_tM = AVCONM+(ORATE_tD-AVCONM)/(1+Ta3), DWIPauxM =
AVCONM*Tpbar3, COMRATE_tM = ORATEaux3M, AINV_tM = AINVM+COMRATE_tM-ORATE_tD, 
EINVauxM = DINVM-AINVM, ORATE_tM = AVCONM+EINVauxM/Ti3+(DWIPauxM-WIPM)/Tw3],(
"Parameters")=[CONSpar = 0, Ta = 8, Tpbar = 2, Ti = 4, Tw = 4, Ta2 = 8, Tpbar2
= 3, Ti2 = 4, Tw2 = 4, Ta3 = 8, Tpbar3 = 4, Ti3 = 4, Tw3 = 4],("AlgEqns")=[],(
"AlgVars")=[],("DynVars")=[AVCON, CONSt, CONSaux, ORATEaux1, COMRATE, ORATE, 
AINV, AVCOND, ORATEaux1D, ORATEaux2D, COMRATED, ORATED, AINVD, AVCONM, 
ORATEaux1M, ORATEaux2M, ORATEaux3M, COMRATEM, ORATEM, AINVM],("DynEqns")=[
AVCON = AVCON_t, CONSt = CONSpar, CONSaux = CONSpar, ORATEaux1 = ORATE, 
COMRATE = COMRATE_t, ORATE = ORATE_t, AINV = AINV_t, AVCOND = AVCON_tD, 
ORATEaux1D = ORATED, ORATEaux2D = ORATEaux1D, COMRATED = COMRATE_tD, ORATED =
ORATE_tD, AINVD = AINV_tD, AVCONM = AVCON_tM, ORATEaux1M = ORATEM, ORATEaux2M
= ORATEaux1M, ORATEaux3M = ORATEaux2M, COMRATEM = COMRATE_tM, ORATEM = 
ORATE_tM, AINVM = AINV_tM]]); getSys := proc () return eval(Syst) end proc; 
end module; end module; Goswami1996 := module () export Sys; Sys := module ()
local aSys, M1, N1, gV1, M1m1, Part2Mgen, Qm, Qp, Mgen1, switchConds, Qpm1, 
Qgen, Tr1, TransformVars1, aPoincSect; export getSys, getSwitchingConditions,
getPoincareSect; aSys := table([("ExplicitAEs")=[beta = 1, g = 9.8, a = 1/(1+
beta), alpha = 1/2*ts-1/2*tn],("Parameters")=[mu = 2, anglePoincare = 3],(
"AlgVars")=[],("DynVars")=[tn, ts, dn, ds],("ODEs")=[`tn'` = dn, `ts'` = ds, 
`dn'` = (mu+2+2*beta*mu+2*beta+beta^2*mu+beta^2)/beta^2/(-mu-2-2*beta*mu-2*
beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)
^2)*(-(1+beta)*beta*ds^2*sin(-ts+tn)+1/a*g*beta*sin(tn))+(1+beta)/beta*cos(2*
alpha)/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*
alpha)^2+beta^2*cos(2*alpha)^2)*((1+beta)*beta*dn^2*sin(-ts+tn)-1/a*((mu+1)*(1
+beta)+1)*g*sin(ts)), `ds'` = (1+beta)/beta*cos(2*alpha)/(-mu-2-2*beta*mu-2*
beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)
^2)*(-(1+beta)*beta*ds^2*sin(-ts+tn)+1/a*g*beta*sin(tn))+1/(-mu-2-2*beta*mu-2*
beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)
^2)*((1+beta)*beta*dn^2*sin(-ts+tn)-1/a*((mu+1)*(1+beta)+1)*g*sin(ts))],("AEs"
)=[]]); M1 := Matrix(2, 2, [[beta^2,-(1+beta)*beta*cos(2*alpha)],[-(1+beta)*
beta*cos(2*alpha),(1+beta)^2*(mu+1)+1]]); N1 := Matrix(2, 2, [[0,-(1+beta)*
beta*ds*sin(-ts+tn)],[(1+beta)*beta*dn*sin(-ts+tn),0]]); gV1 := Vector(2, [g*
beta*sin(tn),-((mu+1)*(1+beta)+1)*g*sin(ts)]); M1m1 := Matrix(2, 2, [[-(mu+2+2
*beta*mu+2*beta+beta^2*mu+beta^2)/beta^2/(-mu-2-2*beta*mu-2*beta-beta^2*mu-
beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2),-(1+beta)/
beta*cos(2*alpha)/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*
beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2)],[-(1+beta)/beta*cos(2*alpha)/(-mu-\
2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^
2*cos(2*alpha)^2),-1/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2
*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2)]]); Part2Mgen := Vector(2, [-(1+
beta)*beta*ds^2*sin(-ts+tn)+1/a*g*beta*sin(tn),(1+beta)*beta*dn^2*sin(-ts+tn)-
1/a*((mu+1)*(1+beta)+1)*g*sin(ts)]); Qm := Matrix(2, 2, [[-beta,-beta+(mu*(1+
beta)^2+2+2*beta)*cos(2*alpha)],[0,-beta]]); Qp := Matrix(2, 2, [[beta*(beta-(
1+beta)*cos(2*alpha)),(1+beta)*(1+beta-beta*cos(2*alpha))+1+mu*(1+beta)^2],[
beta^2,-(1+beta)*beta*cos(2*alpha)]]); Mgen1 := Vector(2, [(mu+2+2*beta*mu+2*
beta+beta^2*mu+beta^2)/beta^2/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*
alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2)*(-(1+beta)*beta*ds^2*sin
(-ts+tn)+1/a*g*beta*sin(tn))+(1+beta)/beta*cos(2*alpha)/(-mu-2-2*beta*mu-2*
beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)
^2)*((1+beta)*beta*dn^2*sin(-ts+tn)-1/a*((mu+1)*(1+beta)+1)*g*sin(ts)),(1+beta
)/beta*cos(2*alpha)/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*
beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2)*(-(1+beta)*beta*ds^2*sin(-ts+tn)+1/
a*g*beta*sin(tn))+1/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*
beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2)*((1+beta)*beta*dn^2*sin(-ts+tn)-1/a
*((mu+1)*(1+beta)+1)*g*sin(ts))]); switchConds := [ts, tn, 2*cos(ts-tn)/(-4*mu
-5+4*cos(ts-tn)^2)*dn+(-2*cos(ts-tn)/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos
(ts-tn))-(-5+2*cos(ts-tn)-4*mu)/(-4*mu-5+4*cos(ts-tn)^2))*ds, 1/(-4*mu-5+4*cos
(ts-tn)^2)*dn+(-1/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos(ts-tn))+(-1+2*cos(
ts-tn))/(-4*mu-5+4*cos(ts-tn)^2))*ds]; Qpm1 := Matrix(2, 2, [[-(1+beta)/beta*
cos(2*alpha)/(-mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*
cos(2*alpha)^2+beta^2*cos(2*alpha)^2),(-2-2*beta+beta*cos(2*alpha)-beta^2+beta
^2*cos(2*alpha)-mu-2*beta*mu-beta^2*mu)/beta^2/(-mu-2-2*beta*mu-2*beta-beta^2*
mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*cos(2*alpha)^2)],[-1/(-
mu-2-2*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+
beta^2*cos(2*alpha)^2),-1/beta*(-beta+cos(2*alpha)+beta*cos(2*alpha))/(-mu-2-2
*beta*mu-2*beta-beta^2*mu-beta^2+cos(2*alpha)^2+2*beta*cos(2*alpha)^2+beta^2*
cos(2*alpha)^2)]]); Qgen := Matrix(2, 2, [[2*cos(ts-tn)/(-4*mu-5+4*cos(ts-tn)^
2),-2*cos(ts-tn)/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos(ts-tn))-(-5+2*cos(
ts-tn)-4*mu)/(-4*mu-5+4*cos(ts-tn)^2)],[1/(-4*mu-5+4*cos(ts-tn)^2),-1/(-4*mu-5
+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos(ts-tn))+(-1+2*cos(ts-tn))/(-4*mu-5+4*cos(ts-
tn)^2)]]); Tr1 := Matrix(4, 4, [[0,1,0,0],[1,0,0,0],[0,0,2*cos(ts-tn)/(-4*mu-5
+4*cos(ts-tn)^2),-2*cos(ts-tn)/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos(ts-tn
))-(-5+2*cos(ts-tn)-4*mu)/(-4*mu-5+4*cos(ts-tn)^2)],[0,0,1/(-4*mu-5+4*cos(ts-
tn)^2),-1/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos(ts-tn))+(-1+2*cos(ts-tn))/
(-4*mu-5+4*cos(ts-tn)^2)]]); TransformVars1 := Vector(4, [ts,tn,2*cos(ts-tn)/(
-4*mu-5+4*cos(ts-tn)^2)*dn+(-2*cos(ts-tn)/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4
)*cos(ts-tn))-(-5+2*cos(ts-tn)-4*mu)/(-4*mu-5+4*cos(ts-tn)^2))*ds,1/(-4*mu-5+4
*cos(ts-tn)^2)*dn+(-1/(-4*mu-5+4*cos(ts-tn)^2)*(-1+(4*mu+4)*cos(ts-tn))+(-1+2*
cos(ts-tn))/(-4*mu-5+4*cos(ts-tn)^2))*ds]); aPoincSect := ts+tn+.3490658504e-1
*anglePoincare; getSys := proc () return eval(aSys) end proc; 
getSwitchingConditions := proc () return eval(switchConds) end proc; 
getPoincareSect := proc () return eval(aPoincSect) end proc; end module; end 
module; John1994 := module () export Sys; Sys := module () local Syst; export
getSys, getFixPoint; Syst := table([("ExplicitAlgEqns")=[CONS_t = CONSpar, 
AVCON_t = AVCON+1/(1+Ta)*(CONS_t-AVCON), DWIPaux = AVCON*Tpbar, EWIPaux = 
DWIPaux-WIP, EINVaux = DINV-AINV, ORATE_t = AVCON_t+EINVaux/Ti+EWIPaux/Tw, 
WIP_t = Tp*(WIP+ORATE_t)/(1+Tp), COMRATE_t = WIP_t/Tp, AINV_t = AINV+COMRATE_t
-CONS_t],("Parameters")=[CONSpar = 4, Ti = 16, Ta = 16, Tw = 16, Tp = 4, Tpbar
= 4, DINV = 50],("AlgEqns")=[],("AlgVars")=[],("DynVars")=[AINV, WIP, AVCON],(
"DynEqns")=[AINV = AINV_t, WIP = WIP_t, AVCON = AVCON_t]]); getSys := proc ()
return eval(Syst) end proc; getFixPoint := proc () local FP; FP := [Ti = 3., 
Ta = 1., Tw = 16., Tp = 3., Tpbar = 3., AINV = 1., DINV = 1., CONSpar = 1., 
AVCON = 1., WIP = 3.]; return FP end proc; end module; end module; Lehman1994
:= module () export Sys; Sys := module () local aSys; export getSys; aSys := 
table([("ExplicitAEs")=[B = 20, beta = 2.4, gamma1 = 12, lambda = .25, x2c = 0
, x10 = 1],("Parameters")=[x20 = 3, Da = .2],("AlgVars")=[],("DynVars")=[x1, 
x2],("ODEs")=[`x1'` = x10-1/lambda*x1+(1/lambda-1)*x1Del-Da*exp(x2/(1+x2/
gamma1))*x1, `x2'` = x20-1/lambda*x2+(1/lambda-1)*x2Del+B*Da*exp(x2/(1+x2/
gamma1))*x1-beta*(x2-x2c)],("DelVars")=[x1Del, x2Del],("AEs")=[]]); getSys :=
proc () return eval(aSys) end proc; end module; end module; Luzyanina2001 := 
module () export Sys; Sys := module () local aSys; export getSys; aSys := 
table([("ExplicitAEs")=[],("Parameters")=[k = .5, alpha = .5],("AlgVars")=[],(
"DynVars")=[x1],("ODEs")=[`x1'` = -(k+alpha^2)*x1Del],("DelVars")=[x1Del],(
"AEs")=[]]); getSys := proc () return eval(aSys) end proc; end module; end 
module; Minimal := module () export Sys; Sys := module () local aSys; export 
getSys; aSys := table([("ExplicitAEs")=[gamma = 1],("Parameters")=[alpha = 1,
beta = 1],("Delays")=[tau[1] = 1],("AlgVars")=[],("DynVars")=[x1, x2],("ODEs")
=[`x1'` = alpha*x1Del-gamma*exp(x2), `x2'` = -beta*x1],("DelVars")=[x1Del, 
x2Del],("AEs")=[]]); getSys := proc () return eval(aSys) end proc; end module;
end module; Olsen1983 := module () export Sys; Sys := module () local aSys; 
export getSys; aSys := table([("ExplicitAEs")=[k2 = 250, k4 = 20, k5 = 5.35, 
k6 = 1/100000, k7 = .1, k8 = .825, A0 = 8],("Parameters")=[k3 = .1, k1 = .3],(
"AlgVars")=[],("DynVars")=[A, B, X, Y],("ODEs")=[`A'` = k7*(A0-A)-k3*A*B*Y, 
`B'` = k8-k1*B*X-k3*A*B*Y, `X'` = k1*B*X-2*k2*X^2+3*k3*A*B*Y-k4*X+k6, `Y'` = 2
*k2*X^2-k5*Y-k3*A*B*Y],("AEs")=[]]); getSys := proc () return eval(aSys) end 
proc; end module; end module; Scott1990 := module () export Sys; Sys := module
() local aSys; export getSys; aSys := table([("ExplicitAEs")=[k = .5500000000e\
-2, delta = .1],("Parameters")=[mu0 = .494710, g = .249179],("AlgVars")=[],(
"DynVars")=[a, b, phi],("ODEs")=[`a'` = mu0*exp(phi*delta)-a*b^2-k*a, `b'` = a
*b^2+k*a-b, `phi'` = b-g*phi],("AEs")=[]]); getSys := proc () return eval(aSys
) end proc; end module; end module; SOFC := module () export Reformer; 
Reformer := module () local Sys; export getSys, getSol; Sys := table([(
"ExplicitAEs")=[r1 = .175e-1, r2 = .185e-1, r3 = .35e-1, r4 = .36e-1, dh_POX =
2*r3-2*r2, dh_NB = 2*r1, rho_POXwand = 7900, rho_NBwand = 7900, C_POXwand = 
477, C_NBwand = 477, Sigma = .567e-7, Epsilon_NBwanda = .9, Epsilon_POXwand =
.9, por = 1, LPOX = .15, A_NBwanda = 2*Pi*r2*LPOX, A_POXwand = 2*Pi*r3*LPOX, 
A_NBwandi = 2*Pi*r1*LPOX, V_POX = LPOX*Pi*(r3^2-r2^2), V_NB = LPOX*Pi*r1^2, R
= 8.31451, p0 = 101300, Mi_H2 = .2e-2, Mi_CO = .28e-1, Mi_C3H8 = .44e-1, 
Mi_CH4 = .16e-1, Mi_H2O = .18e-1, Mi_CO2 = .44e-1, Mi_N2 = .28e-1, Mi_O2 = .32\
e-1, C = r3/r2, B = (C^2-1)/ln(C), D = 1+C^2-B, Nu_Corr = 8, Nu_NBwanda = 
Nu_Corr*D*(1-C)/C/(3/8+7/8*C^2-1/2*B-1/2*C^4/B), Nu_POXwand = Nu_Corr*D*(1-C)/
C/(1/2*B-7/8-3/8*C^2+1/2/B), Nu_NBwandi = 4.2*Nu_Corr, xi_POXLuftH2 = 0, 
xi_POXLuftCO = 0, xi_POXLuftC3H8 = 0, xi_POXLuftCH4 = 0, xi_POXLuftH2O = 0, 
xi_POXLuftCO2 = 0, xi_POXLuftN2 = .79, xi_POXLuftO2 = .21, xi_POXeH2 = 
ni_POXeH2/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+
ni_POXeN2+ni_POXeO2), xi_POXeCO = ni_POXeCO/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8+
ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+ni_POXeN2+ni_POXeO2), xi_POXeC3H8 = 
ni_POXeC3H8/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+
ni_POXeN2+ni_POXeO2), xi_POXeCH4 = ni_POXeCH4/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8
+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+ni_POXeN2+ni_POXeO2), xi_POXeH2O = 
ni_POXeH2O/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+
ni_POXeN2+ni_POXeO2), xi_POXeCO2 = ni_POXeCO2/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8
+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+ni_POXeN2+ni_POXeO2), xi_POXeN2 = ni_POXeN2/
(ni_POXeH2+ni_POXeCO+ni_POXeC3H8+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+ni_POXeN2+
ni_POXeO2), xi_POXeO2 = ni_POXeO2/(ni_POXeH2+ni_POXeCO+ni_POXeC3H8+ni_POXeCH4+
ni_POXeH2O+ni_POXeCO2+ni_POXeN2+ni_POXeO2), xi_NBeH2 = xi_anaH2, xi_NBeCO = 
xi_anaCO, xi_NBeC3H8 = xi_anaC3H8, xi_NBeCH4 = xi_anaCH4, xi_NBeH2O = 
xi_anaH2O, xi_NBeCO2 = xi_anaCO2, xi_NBeN2 = xi_anaN2, xi_NBeO2 = xi_anaO2, 
xi_POXLuftH2 = 0, xi_POXLuftCO = 0, xi_POXLuftC3H8 = 0, xi_POXLuftCH4 = 0, 
xi_POXLuftH2O = 0, xi_POXLuftCO2 = 0, xi_POXLuftN2 = .79, xi_POXLuftO2 = .21,
xi_NBe2H2 = xi_anaH2, xi_NBe2CO = xi_anaCO, xi_NBe2C3H8 = xi_anaC3H8, 
xi_NBe2CH4 = xi_anaCH4, xi_NBe2H2O = xi_anaH2O, xi_NBe2CO2 = xi_anaCO2, 
xi_NBe2N2 = xi_anaN2, xi_NBe2O2 = xi_anaO2, ni_POXmH2 = ni_POXeH2+ni_POXLuftH2
, ni_POXmCO = ni_POXeCO+ni_POXLuftCO, ni_POXmC3H8 = ni_POXeC3H8+ni_POXLuftC3H8
, ni_POXmCH4 = ni_POXeCH4+ni_POXLuftCH4, ni_POXmH2O = ni_POXeH2O+ni_POXLuftH2O
, ni_POXmCO2 = ni_POXeCO2+ni_POXLuftCO2, ni_POXmN2 = ni_POXeN2+ni_POXLuftN2, 
ni_POXmO2 = ni_POXeO2+ni_POXLuftO2, ni_NBmH2 = ni_NBeH2+ni_NBe2H2, ni_NBmCO =
ni_NBeCO+ni_NBe2CO, ni_NBmC3H8 = ni_NBeC3H8+ni_NBe2C3H8, ni_NBmCH4 = ni_NBeCH4
+ni_NBe2CH4, ni_NBmH2O = ni_NBeH2O+ni_NBe2H2O, ni_NBmCO2 = ni_NBeCO2+
ni_NBe2CO2, ni_NBmN2 = ni_NBeN2+ni_NBe2N2, ni_NBmO2 = ni_NBeO2+ni_NBe2O2, 
xi_POXmH2 = ni_POXmH2/(ni_POXmH2+ni_POXmCO+ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+
ni_POXmCO2+ni_POXmN2+ni_POXmO2), xi_POXmCO = ni_POXmCO/(ni_POXmH2+ni_POXmCO+
ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+ni_POXmCO2+ni_POXmN2+ni_POXmO2), xi_POXmC3H8
= ni_POXmC3H8/(ni_POXmH2+ni_POXmCO+ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+
ni_POXmCO2+ni_POXmN2+ni_POXmO2), xi_POXmCH4 = ni_POXmCH4/(ni_POXmH2+ni_POXmCO+
ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+ni_POXmCO2+ni_POXmN2+ni_POXmO2), xi_POXmH2O
= ni_POXmH2O/(ni_POXmH2+ni_POXmCO+ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+ni_POXmCO2
+ni_POXmN2+ni_POXmO2), xi_POXmCO2 = ni_POXmCO2/(ni_POXmH2+ni_POXmCO+
ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+ni_POXmCO2+ni_POXmN2+ni_POXmO2), xi_POXmN2 =
ni_POXmN2/(ni_POXmH2+ni_POXmCO+ni_POXmC3H8+ni_POXmCH4+ni_POXmH2O+ni_POXmCO2+
ni_POXmN2+ni_POXmO2), xi_POXmO2 = ni_POXmO2/(ni_POXmH2+ni_POXmCO+ni_POXmC3H8+
ni_POXmCH4+ni_POXmH2O+ni_POXmCO2+ni_POXmN2+ni_POXmO2), xi_POXaC3H8 = 0, 
xi_POXaO2 = 0, m_POXm = Mi_H2*ni_POXmH2+Mi_CO*ni_POXmCO+Mi_C3H8*ni_POXmC3H8+
Mi_CH4*ni_POXmCH4+Mi_H2O*ni_POXmH2O+Mi_CO2*ni_POXmCO2+Mi_N2*ni_POXmN2+Mi_O2*
ni_POXmO2, n_POXa = m_POXm/(Mi_H2*xi_POXaH2+Mi_CO*xi_POXaCO+Mi_C3H8*
xi_POXaC3H8+Mi_CH4*xi_POXaCH4+Mi_H2O*xi_POXaH2O+Mi_CO2*xi_POXaCO2+Mi_N2*
xi_POXaN2+Mi_O2*xi_POXaO2), ni_POXaH2 = xi_POXaH2*n_POXa, ni_POXaCO = 
xi_POXaCO*n_POXa, ni_POXaCH4 = xi_POXaCH4*n_POXa, ni_POXaC3H8 = xi_POXaC3H8*
n_POXa, ni_POXaCO2 = xi_POXaCO2*n_POXa, ni_POXaH2O = xi_POXaH2O*n_POXa, 
ni_POXaN2 = xi_POXaN2*n_POXa, ni_POXaO2 = xi_POXaO2*n_POXa, ni_NBrH2 = 
ni_NBmH2, ni_NBrCO = ni_NBmCO, ni_NBrC3H8 = ni_NBmC3H8, ni_NBrCH4 = ni_NBmCH4,
ni_NBrO2 = .5*ni_NBrH2+.5*ni_NBrCO+2*ni_NBrCH4+5*ni_NBrC3H8, ni_NBrCO2 = 
ni_NBrCO+ni_NBrCH4+3*ni_NBrC3H8, ni_NBrH2O = ni_NBrH2+2*ni_NBrCH4+4*ni_NBrC3H8
, ni_NBaH2 = ni_NBmH2-ni_NBrH2, ni_NBaCO = ni_NBmCO-ni_NBrCO, ni_NBaCH4 = 
ni_NBmCH4-ni_NBrCH4, ni_NBaC3H8 = ni_NBmC3H8-ni_NBrC3H8, ni_NBaCO2 = ni_NBmCO2
+ni_NBrCO2, ni_NBaH2O = ni_NBmH2O+ni_NBrH2O, ni_NBaN2 = ni_NBmN2+ni_NBLuftN2,
ni_NBaO2 = ni_NBmO2+ni_NBLuftO2-ni_NBrO2, xi_NBaH2 = ni_NBaH2/(ni_NBaH2+
ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+ni_NBaN2+ni_NBaO2), xi_NBaCO
= ni_NBaCO/(ni_NBaH2+ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+
ni_NBaN2+ni_NBaO2), xi_NBaC3H8 = ni_NBaC3H8/(ni_NBaH2+ni_NBaCO+ni_NBaC3H8+
ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+ni_NBaN2+ni_NBaO2), xi_NBaCH4 = ni_NBaCH4/(
ni_NBaH2+ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+ni_NBaN2+ni_NBaO2),
xi_NBaH2O = ni_NBaH2O/(ni_NBaH2+ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+
ni_NBaCO2+ni_NBaN2+ni_NBaO2), xi_NBaCO2 = ni_NBaCO2/(ni_NBaH2+ni_NBaCO+
ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+ni_NBaN2+ni_NBaO2), xi_NBaN2 = 
ni_NBaN2/(ni_NBaH2+ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+ni_NBaN2+
ni_NBaO2), xi_NBaO2 = ni_NBaO2/(ni_NBaH2+ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+
ni_NBaH2O+ni_NBaCO2+ni_NBaN2+ni_NBaO2), p_POXa = 1/2*p_POXe+1/2*p_POXLuft, 
p_NBa = 1/3*p_NBe+1/3*p_NBLuft+1/3*p_NBe2, h_POXeH2 = (-40783.23210/T_POXe^2-\
800.9186040/T_POXe*ln(T_POXe)+8.214702010-.6348572285e-2*T_POXe+.5845350253e-5
*T_POXe^2-.3007150675e-8*T_POXe^3+.6736186980e-12*T_POXe^4+2682.484665/T_POXe)
*T_POXe*R, h_POXeCO = (-14890.45326/T_POXe^2-292.2285939/T_POXe*ln(T_POXe)+5.7\
24527170-.4088117515e-2*T_POXe+.4856344897e-5*T_POXe^2-.2719365755e-8*T_POXe^3
+.6055883654e-12*T_POXe^4-13031.31878/T_POXe)*T_POXe*R, h_POXeC3H8 = (243314.4\
337/T_POXe^2+4656.270810/T_POXe*ln(T_POXe)-29.39466091+.5944763725e-1*T_POXe-.\
4587694230e-4*T_POXe^2+.2203705978e-7*T_POXe^3-.4685975988e-11*T_POXe^4-35403.\
35270/T_POXe)*T_POXe*R, h_POXeCH4 = (176685.0998/T_POXe^2+2786.181020/T_POXe*
ln(T_POXe)-12.02577850+.1958809645e-1*T_POXe-.1206351477e-4*T_POXe^2+.\
5067132608e-8*T_POXe^3-.9953410980e-12*T_POXe^4-23313.14360/T_POXe)*T_POXe*R,
h_POXeH2O = (39479.60830/T_POXe^2+575.5731020/T_POXe*ln(T_POXe)+.9317826530+.\
3611356430e-2*T_POXe-.2447519123e-5*T_POXe^2+.1238760872e-8*T_POXe^3-.\
2673866492e-12*T_POXe^4-33039.74310/T_POXe)*T_POXe*R, h_POXeCO2 = (-49436.5054\
0/T_POXe^2-626.4116010/T_POXe*ln(T_POXe)+5.301725240+.1251906908e-2*T_POXe-.\
7091029093e-7*T_POXe^2-.1922497195e-9*T_POXe^3+.5699355602e-13*T_POXe^4-45281.\
98460/T_POXe)*T_POXe*R, h_POXeN2 = (-22103.71497/T_POXe^2-381.8461820/T_POXe*
ln(T_POXe)+6.082738360-.4265457205e-2*T_POXe+.4615487297e-5*T_POXe^2-.\
2406448405e-8*T_POXe^3+.5039411618e-12*T_POXe^4+710.8460860/T_POXe)*T_POXe*R,
h_POXeO2 = (34255.63420/T_POXe^2+484.7000970/T_POXe*ln(T_POXe)+1.119010961+.\
2146944620e-2*T_POXe-.2278766840e-6*T_POXe^2-.5058431750e-9*T_POXe^3+.\
2078080036e-12*T_POXe^4-3391.454870/T_POXe)*T_POXe*R, h_POXLuftH2 = (-40783.23\
210/T_POXLuft^2-800.9186040/T_POXLuft*ln(T_POXLuft)+8.214702010-.6348572285e-2
*T_POXLuft+.5845350253e-5*T_POXLuft^2-.3007150675e-8*T_POXLuft^3+.6736186980e-\
12*T_POXLuft^4+2682.484665/T_POXLuft)*T_POXLuft*R, h_POXLuftCO = (-14890.45326
/T_POXLuft^2-292.2285939/T_POXLuft*ln(T_POXLuft)+5.724527170-.4088117515e-2*
T_POXLuft+.4856344897e-5*T_POXLuft^2-.2719365755e-8*T_POXLuft^3+.6055883654e-\
12*T_POXLuft^4-13031.31878/T_POXLuft)*T_POXLuft*R, h_POXLuftC3H8 = (243314.433\
7/T_POXLuft^2+4656.270810/T_POXLuft*ln(T_POXLuft)-29.39466091+.5944763725e-1*
T_POXLuft-.4587694230e-4*T_POXLuft^2+.2203705978e-7*T_POXLuft^3-.4685975988e-\
11*T_POXLuft^4-35403.35270/T_POXLuft)*T_POXLuft*R, h_POXLuftCH4 = (176685.0998
/T_POXLuft^2+2786.181020/T_POXLuft*ln(T_POXLuft)-12.02577850+.1958809645e-1*
T_POXLuft-.1206351477e-4*T_POXLuft^2+.5067132608e-8*T_POXLuft^3-.9953410980e-\
12*T_POXLuft^4-23313.14360/T_POXLuft)*T_POXLuft*R, h_POXLuftH2O = (39479.60830
/T_POXLuft^2+575.5731020/T_POXLuft*ln(T_POXLuft)+.9317826530+.3611356430e-2*
T_POXLuft-.2447519123e-5*T_POXLuft^2+.1238760872e-8*T_POXLuft^3-.2673866492e-\
12*T_POXLuft^4-33039.74310/T_POXLuft)*T_POXLuft*R, h_POXLuftCO2 = (-49436.5054\
0/T_POXLuft^2-626.4116010/T_POXLuft*ln(T_POXLuft)+5.301725240+.1251906908e-2*
T_POXLuft-.7091029093e-7*T_POXLuft^2-.1922497195e-9*T_POXLuft^3+.5699355602e-\
13*T_POXLuft^4-45281.98460/T_POXLuft)*T_POXLuft*R, h_POXLuftN2 = (-22103.71497
/T_POXLuft^2-381.8461820/T_POXLuft*ln(T_POXLuft)+6.082738360-.4265457205e-2*
T_POXLuft+.4615487297e-5*T_POXLuft^2-.2406448405e-8*T_POXLuft^3+.5039411618e-\
12*T_POXLuft^4+710.8460860/T_POXLuft)*T_POXLuft*R, h_POXLuftO2 = (34255.63420/
T_POXLuft^2+484.7000970/T_POXLuft*ln(T_POXLuft)+1.119010961+.2146944620e-2*
T_POXLuft-.2278766840e-6*T_POXLuft^2-.5058431750e-9*T_POXLuft^3+.2078080036e-\
12*T_POXLuft^4-3391.454870/T_POXLuft)*T_POXLuft*R, h_POXaH2 = (-560812.8010/
T_POXa^2-837.1504740/T_POXa*ln(T_POXa)+2.975364532+.6261245620e-3*T_POXa-.\
1246905397e-6*T_POXa^2+.1484156300e-10*T_POXa^3-.7213988200e-15*T_POXa^4+5339.\
824410/T_POXa)*T_POXa*R, h_POXaCO = (-461919.7250/T_POXa^2-1944.704863/T_POXa*
ln(T_POXa)+5.916714180-.2832141415e-3*T_POXa+.4662715133e-7*T_POXa^2-.\
4469200902e-11*T_POXa^3+.1924187114e-15*T_POXa^4-2466.261084/T_POXa)*T_POXa*R,
h_POXaC3H8 = (-6420731.680/T_POXa^2-26597.91134/T_POXa*ln(T_POXa)+45.34356840-
.2510331960e-2*T_POXa+.3157072313e-6*T_POXa^2-.2393851308e-10*T_POXa^3+.\
8019345760e-15*T_POXa^4+145558.2459/T_POXa)*T_POXa*R, h_POXaCH4 = (-3730042.76\
0/T_POXa^2-13835.01485/T_POXa*ln(T_POXa)+20.49107091-.9809873795e-3*T_POXa+.\
1575771013e-6*T_POXa^2-.9322036725e-11*T_POXa^3+.3247474414e-15*T_POXa^4+75320\
.66910/T_POXa)*T_POXa*R, h_POXaH2O = (-1034972.096/T_POXa^2-2412.698562/T_POXa
*ln(T_POXa)+4.646110780+.1145999154e-2*T_POXa-.2278943493e-6*T_POXa^2+.\
2356617232e-10*T_POXa^3-.9644761060e-15*T_POXa^4-13842.86509/T_POXa)*T_POXa*R,
h_POXaCO2 = (-117696.2419/T_POXa^2-1788.791477/T_POXa*ln(T_POXa)+8.291523190-.\
4611578390e-4*T_POXa+.1621225627e-8*T_POXa^2-.4727633280e-12*T_POXa^3+.\
1266007318e-15*T_POXa^4-39083.50590/T_POXa)*T_POXa*R, h_POXaN2 = (-587712.4060
/T_POXa^2-2239.249073/T_POXa*ln(T_POXa)+6.066949220-.3069842750e-3*T_POXa+.\
4972688930e-7*T_POXa^2-.4807763712e-11*T_POXa^3+.2123908772e-15*T_POXa^4+12832\
.10415/T_POXa)*T_POXa*R, h_POXaO2 = (1037939.022/T_POXa^2+2344.830282/T_POXa*
ln(T_POXa)+1.819732036+.6339237910e-3*T_POXa-.7293559960e-7*T_POXa^2+.\
5134298930e-11*T_POXa^3-.1638693410e-15*T_POXa^4-16890.10929/T_POXa)*T_POXa*R,
Cp_POXaH2 = R*(560812.8010/T_POXa^2-837.1504740/T_POXa+2.975364532+.1252249124\
e-2*T_POXa-.3740716190e-6*T_POXa^2+.5936625200e-10*T_POXa^3-.3606994100e-14*
T_POXa^4), Cp_POXaCO = R*(461919.7250/T_POXa^2-1944.704863/T_POXa+5.916714180-
.5664282830e-3*T_POXa+.1398814540e-6*T_POXa^2-.1787680361e-10*T_POXa^3+.\
9620935570e-15*T_POXa^4), Cp_POXaC3H8 = R*(6420731.680/T_POXa^2-26597.91134/
T_POXa+45.34356840-.5020663920e-2*T_POXa+.9471216940e-6*T_POXa^2-.9575405230e-\
10*T_POXa^3+.4009672880e-14*T_POXa^4), Cp_POXaCH4 = R*(3730042.760/T_POXa^2-\
13835.01485/T_POXa+20.49107091-.1961974759e-2*T_POXa+.4727313040e-6*T_POXa^2-.\
3728814690e-10*T_POXa^3+.1623737207e-14*T_POXa^4), Cp_POXaH2O = R*(1034972.096
/T_POXa^2-2412.698562/T_POXa+4.646110780+.2291998307e-2*T_POXa-.6836830480e-6*
T_POXa^2+.9426468930e-10*T_POXa^3-.4822380530e-14*T_POXa^4), Cp_POXaCO2 = R*(
117696.2419/T_POXa^2-1788.791477/T_POXa+8.291523190-.9223156780e-4*T_POXa+.\
4863676880e-8*T_POXa^2-.1891053312e-11*T_POXa^3+.6330036590e-15*T_POXa^4), 
Cp_POXaN2 = R*(587712.4060/T_POXa^2-2239.249073/T_POXa+6.066949220-.6139685500\
e-3*T_POXa+.1491806679e-6*T_POXa^2-.1923105485e-10*T_POXa^3+.1061954386e-14*
T_POXa^4), Cp_POXaO2 = R*(-1037939.022/T_POXa^2+2344.830282/T_POXa+1.819732036
+.1267847582e-2*T_POXa-.2188067988e-6*T_POXa^2+.2053719572e-10*T_POXa^3-.\
8193467050e-15*T_POXa^4), Cp_POXa = xi_POXaH2*Cp_POXaH2+xi_POXaCO*Cp_POXaCO+
xi_POXaC3H8*Cp_POXaC3H8+xi_POXaCH4*Cp_POXaCH4+xi_POXaH2O*Cp_POXaH2O+xi_POXaCO2
*Cp_POXaCO2+xi_POXaN2*Cp_POXaN2+xi_POXaO2*Cp_POXaO2, s_POXaH2 = R*(-280406.400\
5/T_POXa^2+837.1504740/T_POXa+2.975364532*ln(T_POXa)+.1252249124e-2*T_POXa-.\
1870358095e-6*T_POXa^2+.1978875067e-10*T_POXa^3-.9017485250e-15*T_POXa^4-2.202\
774769), s_POXaCO = R*(-230959.8625/T_POXa^2+1944.704863/T_POXa+5.916714180*ln
(T_POXa)-.5664282830e-3*T_POXa+.6994072700e-7*T_POXa^2-.5958934537e-11*T_POXa^
3+.2405233892e-15*T_POXa^4-13.87413108), s_POXaC3H8 = R*(-3210365.840/T_POXa^2
+26597.91134/T_POXa+45.34356840*ln(T_POXa)-.5020663920e-2*T_POXa+.4735608470e-\
6*T_POXa^2-.3191801743e-10*T_POXa^3+.1002418220e-14*T_POXa^4-281.8374734), 
s_POXaCH4 = R*(-1865021.380/T_POXa^2+13835.01485/T_POXa+20.49107091*ln(T_POXa)
-.1961974759e-2*T_POXa+.2363656520e-6*T_POXa^2-.1242938230e-10*T_POXa^3+.\
4059343018e-15*T_POXa^4-121.9124889), s_POXaH2O = R*(-517486.0480/T_POXa^2+
2412.698562/T_POXa+4.646110780*ln(T_POXa)+.2291998307e-2*T_POXa-.3418415240e-6
*T_POXa^2+.3142156310e-10*T_POXa^3-.1205595132e-14*T_POXa^4-7.978148510), 
s_POXaCO2 = R*(-58848.12095/T_POXa^2+1788.791477/T_POXa+8.291523190*ln(T_POXa)
-.9223156780e-4*T_POXa+.2431838440e-8*T_POXa^2-.6303511040e-12*T_POXa^3+.\
1582509148e-15*T_POXa^4-26.52669281), s_POXaN2 = R*(-293856.2030/T_POXa^2+2239\
.249073/T_POXa+6.066949220*ln(T_POXa)-.6139685500e-3*T_POXa+.7459033395e-7*
T_POXa^2-.6410351617e-11*T_POXa^3+.2654885965e-15*T_POXa^4-1.586639599), 
s_POXaO2 = R*(518969.5110/T_POXa^2-2344.830282/T_POXa+1.819732036*ln(T_POXa)+.\
1267847582e-2*T_POXa-.1094033994e-6*T_POXa^2+.6845731907e-11*T_POXa^3-.\
2048366762e-15*T_POXa^4+17.38716506), G_H2 = h_POXaH2-T_POXa*s_POXaH2, G_CO =
h_POXaCO-T_POXa*s_POXaCO, G_CH4 = h_POXaCH4-T_POXa*s_POXaCH4, G_H2O = 
h_POXaH2O-T_POXa*s_POXaH2O, G_CO2 = h_POXaCO2-T_POXa*s_POXaCO2, h_NBeH2 = (-\
560812.8010/T_NBe^2-837.1504740/T_NBe*ln(T_NBe)+2.975364532+.6261245620e-3*
T_NBe-.1246905397e-6*T_NBe^2+.1484156300e-10*T_NBe^3-.7213988200e-15*T_NBe^4+
5339.824410/T_NBe)*T_NBe*R, h_NBeCO = (-461919.7250/T_NBe^2-1944.704863/T_NBe*
ln(T_NBe)+5.916714180-.2832141415e-3*T_NBe+.4662715133e-7*T_NBe^2-.4469200902e\
-11*T_NBe^3+.1924187114e-15*T_NBe^4-2466.261084/T_NBe)*T_NBe*R, h_NBeC3H8 = (-\
6420731.680/T_NBe^2-26597.91134/T_NBe*ln(T_NBe)+45.34356840-.2510331960e-2*
T_NBe+.3157072313e-6*T_NBe^2-.2393851308e-10*T_NBe^3+.8019345760e-15*T_NBe^4+
145558.2459/T_NBe)*T_NBe*R, h_NBeCH4 = (-3730042.760/T_NBe^2-13835.01485/T_NBe
*ln(T_NBe)+20.49107091-.9809873795e-3*T_NBe+.1575771013e-6*T_NBe^2-.9322036725\
e-11*T_NBe^3+.3247474414e-15*T_NBe^4+75320.66910/T_NBe)*T_NBe*R, h_NBeH2O = (-\
1034972.096/T_NBe^2-2412.698562/T_NBe*ln(T_NBe)+4.646110780+.1145999154e-2*
T_NBe-.2278943493e-6*T_NBe^2+.2356617232e-10*T_NBe^3-.9644761060e-15*T_NBe^4-\
13842.86509/T_NBe)*T_NBe*R, h_NBeCO2 = (-117696.2419/T_NBe^2-1788.791477/T_NBe
*ln(T_NBe)+8.291523190-.4611578390e-4*T_NBe+.1621225627e-8*T_NBe^2-.4727633280\
e-12*T_NBe^3+.1266007318e-15*T_NBe^4-39083.50590/T_NBe)*T_NBe*R, h_NBeN2 = (-\
587712.4060/T_NBe^2-2239.249073/T_NBe*ln(T_NBe)+6.066949220-.3069842750e-3*
T_NBe+.4972688930e-7*T_NBe^2-.4807763712e-11*T_NBe^3+.2123908772e-15*T_NBe^4+
12832.10415/T_NBe)*T_NBe*R, h_NBeO2 = (1037939.022/T_NBe^2+2344.830282/T_NBe*
ln(T_NBe)+1.819732036+.6339237910e-3*T_NBe-.7293559960e-7*T_NBe^2+.5134298930e\
-11*T_NBe^3-.1638693410e-15*T_NBe^4-16890.10929/T_NBe)*T_NBe*R, h_NBe2H2 = (-\
40783.23210/T_NBe2^2-800.9186040/T_NBe2*ln(T_NBe2)+8.214702010-.6348572285e-2*
T_NBe2+.5845350253e-5*T_NBe2^2-.3007150675e-8*T_NBe2^3+.6736186980e-12*T_NBe2^
4+2682.484665/T_NBe2)*T_NBe2*R, h_NBe2CO = (-14890.45326/T_NBe2^2-292.2285939/
T_NBe2*ln(T_NBe2)+5.724527170-.4088117515e-2*T_NBe2+.4856344897e-5*T_NBe2^2-.\
2719365755e-8*T_NBe2^3+.6055883654e-12*T_NBe2^4-13031.31878/T_NBe2)*T_NBe2*R,
h_NBe2C3H8 = (243314.4337/T_NBe2^2+4656.270810/T_NBe2*ln(T_NBe2)-29.39466091+.\
5944763725e-1*T_NBe2-.4587694230e-4*T_NBe2^2+.2203705978e-7*T_NBe2^3-.\
4685975988e-11*T_NBe2^4-35403.35270/T_NBe2)*T_NBe2*R, h_NBe2CH4 = (176685.0998
/T_NBe2^2+2786.181020/T_NBe2*ln(T_NBe2)-12.02577850+.1958809645e-1*T_NBe2-.\
1206351477e-4*T_NBe2^2+.5067132608e-8*T_NBe2^3-.9953410980e-12*T_NBe2^4-23313.\
14360/T_NBe2)*T_NBe2*R, h_NBe2H2O = (39479.60830/T_NBe2^2+575.5731020/T_NBe2*
ln(T_NBe2)+.9317826530+.3611356430e-2*T_NBe2-.2447519123e-5*T_NBe2^2+.\
1238760872e-8*T_NBe2^3-.2673866492e-12*T_NBe2^4-33039.74310/T_NBe2)*T_NBe2*R,
h_NBe2CO2 = (-49436.50540/T_NBe2^2-626.4116010/T_NBe2*ln(T_NBe2)+5.301725240+.\
1251906908e-2*T_NBe2-.7091029093e-7*T_NBe2^2-.1922497195e-9*T_NBe2^3+.\
5699355602e-13*T_NBe2^4-45281.98460/T_NBe2)*T_NBe2*R, h_NBe2N2 = (-22103.71497
/T_NBe2^2-381.8461820/T_NBe2*ln(T_NBe2)+6.082738360-.4265457205e-2*T_NBe2+.\
4615487297e-5*T_NBe2^2-.2406448405e-8*T_NBe2^3+.5039411618e-12*T_NBe2^4+710.84\
60860/T_NBe2)*T_NBe2*R, h_NBe2O2 = (34255.63420/T_NBe2^2+484.7000970/T_NBe2*ln
(T_NBe2)+1.119010961+.2146944620e-2*T_NBe2-.2278766840e-6*T_NBe2^2-.5058431750\
e-9*T_NBe2^3+.2078080036e-12*T_NBe2^4-3391.454870/T_NBe2)*T_NBe2*R, h_NBLuftH2
= (-40783.23210/T_NBLuft^2-800.9186040/T_NBLuft*ln(T_NBLuft)+8.214702010-.\
6348572285e-2*T_NBLuft+.5845350253e-5*T_NBLuft^2-.3007150675e-8*T_NBLuft^3+.\
6736186980e-12*T_NBLuft^4+2682.484665/T_NBLuft)*T_NBLuft*R, h_NBLuftCO = (-\
14890.45326/T_NBLuft^2-292.2285939/T_NBLuft*ln(T_NBLuft)+5.724527170-.\
4088117515e-2*T_NBLuft+.4856344897e-5*T_NBLuft^2-.2719365755e-8*T_NBLuft^3+.\
6055883654e-12*T_NBLuft^4-13031.31878/T_NBLuft)*T_NBLuft*R, h_NBLuftC3H8 = (
243314.4337/T_NBLuft^2+4656.270810/T_NBLuft*ln(T_NBLuft)-29.39466091+.\
5944763725e-1*T_NBLuft-.4587694230e-4*T_NBLuft^2+.2203705978e-7*T_NBLuft^3-.\
4685975988e-11*T_NBLuft^4-35403.35270/T_NBLuft)*T_NBLuft*R, h_NBLuftCH4 = (
176685.0998/T_NBLuft^2+2786.181020/T_NBLuft*ln(T_NBLuft)-12.02577850+.\
1958809645e-1*T_NBLuft-.1206351477e-4*T_NBLuft^2+.5067132608e-8*T_NBLuft^3-.\
9953410980e-12*T_NBLuft^4-23313.14360/T_NBLuft)*T_NBLuft*R, h_NBLuftH2O = (
39479.60830/T_NBLuft^2+575.5731020/T_NBLuft*ln(T_NBLuft)+.9317826530+.\
3611356430e-2*T_NBLuft-.2447519123e-5*T_NBLuft^2+.1238760872e-8*T_NBLuft^3-.\
2673866492e-12*T_NBLuft^4-33039.74310/T_NBLuft)*T_NBLuft*R, h_NBLuftCO2 = (-\
49436.50540/T_NBLuft^2-626.4116010/T_NBLuft*ln(T_NBLuft)+5.301725240+.\
1251906908e-2*T_NBLuft-.7091029093e-7*T_NBLuft^2-.1922497195e-9*T_NBLuft^3+.\
5699355602e-13*T_NBLuft^4-45281.98460/T_NBLuft)*T_NBLuft*R, h_NBLuftN2 = (-\
22103.71497/T_NBLuft^2-381.8461820/T_NBLuft*ln(T_NBLuft)+6.082738360-.\
4265457205e-2*T_NBLuft+.4615487297e-5*T_NBLuft^2-.2406448405e-8*T_NBLuft^3+.\
5039411618e-12*T_NBLuft^4+710.8460860/T_NBLuft)*T_NBLuft*R, h_NBLuftO2 = (
34255.63420/T_NBLuft^2+484.7000970/T_NBLuft*ln(T_NBLuft)+1.119010961+.\
2146944620e-2*T_NBLuft-.2278766840e-6*T_NBLuft^2-.5058431750e-9*T_NBLuft^3+.\
2078080036e-12*T_NBLuft^4-3391.454870/T_NBLuft)*T_NBLuft*R, h_NBaH2 = (-560812\
.8010/T_NBa^2-837.1504740/T_NBa*ln(T_NBa)+2.975364532+.6261245620e-3*T_NBa-.\
1246905397e-6*T_NBa^2+.1484156300e-10*T_NBa^3-.7213988200e-15*T_NBa^4+5339.824\
410/T_NBa)*T_NBa*R, h_NBaCO = (-461919.7250/T_NBa^2-1944.704863/T_NBa*ln(T_NBa
)+5.916714180-.2832141415e-3*T_NBa+.4662715133e-7*T_NBa^2-.4469200902e-11*
T_NBa^3+.1924187114e-15*T_NBa^4-2466.261084/T_NBa)*T_NBa*R, h_NBaC3H8 = (-\
6420731.680/T_NBa^2-26597.91134/T_NBa*ln(T_NBa)+45.34356840-.2510331960e-2*
T_NBa+.3157072313e-6*T_NBa^2-.2393851308e-10*T_NBa^3+.8019345760e-15*T_NBa^4+
145558.2459/T_NBa)*T_NBa*R, h_NBaCH4 = (-3730042.760/T_NBa^2-13835.01485/T_NBa
*ln(T_NBa)+20.49107091-.9809873795e-3*T_NBa+.1575771013e-6*T_NBa^2-.9322036725\
e-11*T_NBa^3+.3247474414e-15*T_NBa^4+75320.66910/T_NBa)*T_NBa*R, h_NBaH2O = (-\
1034972.096/T_NBa^2-2412.698562/T_NBa*ln(T_NBa)+4.646110780+.1145999154e-2*
T_NBa-.2278943493e-6*T_NBa^2+.2356617232e-10*T_NBa^3-.9644761060e-15*T_NBa^4-\
13842.86509/T_NBa)*T_NBa*R, h_NBaCO2 = (-117696.2419/T_NBa^2-1788.791477/T_NBa
*ln(T_NBa)+8.291523190-.4611578390e-4*T_NBa+.1621225627e-8*T_NBa^2-.4727633280\
e-12*T_NBa^3+.1266007318e-15*T_NBa^4-39083.50590/T_NBa)*T_NBa*R, h_NBaN2 = (-\
587712.4060/T_NBa^2-2239.249073/T_NBa*ln(T_NBa)+6.066949220-.3069842750e-3*
T_NBa+.4972688930e-7*T_NBa^2-.4807763712e-11*T_NBa^3+.2123908772e-15*T_NBa^4+
12832.10415/T_NBa)*T_NBa*R, h_NBaO2 = (1037939.022/T_NBa^2+2344.830282/T_NBa*
ln(T_NBa)+1.819732036+.6339237910e-3*T_NBa-.7293559960e-7*T_NBa^2+.5134298930e\
-11*T_NBa^3-.1638693410e-15*T_NBa^4-16890.10929/T_NBa)*T_NBa*R, Cp_NBaH2 = R*(
560812.8010/T_NBa^2-837.1504740/T_NBa+2.975364532+.1252249124e-2*T_NBa-.\
3740716190e-6*T_NBa^2+.5936625200e-10*T_NBa^3-.3606994100e-14*T_NBa^4), 
Cp_NBaCO = R*(461919.7250/T_NBa^2-1944.704863/T_NBa+5.916714180-.5664282830e-3
*T_NBa+.1398814540e-6*T_NBa^2-.1787680361e-10*T_NBa^3+.9620935570e-15*T_NBa^4)
, Cp_NBaC3H8 = R*(6420731.680/T_NBa^2-26597.91134/T_NBa+45.34356840-.\
5020663920e-2*T_NBa+.9471216940e-6*T_NBa^2-.9575405230e-10*T_NBa^3+.4009672880\
e-14*T_NBa^4), Cp_NBaCH4 = R*(3730042.760/T_NBa^2-13835.01485/T_NBa+20.4910709\
1-.1961974759e-2*T_NBa+.4727313040e-6*T_NBa^2-.3728814690e-10*T_NBa^3+.\
1623737207e-14*T_NBa^4), Cp_NBaH2O = R*(1034972.096/T_NBa^2-2412.698562/T_NBa+
4.646110780+.2291998307e-2*T_NBa-.6836830480e-6*T_NBa^2+.9426468930e-10*T_NBa^
3-.4822380530e-14*T_NBa^4), Cp_NBaCO2 = R*(117696.2419/T_NBa^2-1788.791477/
T_NBa+8.291523190-.9223156780e-4*T_NBa+.4863676880e-8*T_NBa^2-.1891053312e-11*
T_NBa^3+.6330036590e-15*T_NBa^4), Cp_NBaN2 = R*(587712.4060/T_NBa^2-2239.24907\
3/T_NBa+6.066949220-.6139685500e-3*T_NBa+.1491806679e-6*T_NBa^2-.1923105485e-\
10*T_NBa^3+.1061954386e-14*T_NBa^4), Cp_NBaO2 = R*(-1037939.022/T_NBa^2+2344.8\
30282/T_NBa+1.819732036+.1267847582e-2*T_NBa-.2188067988e-6*T_NBa^2+.\
2053719572e-10*T_NBa^3-.8193467050e-15*T_NBa^4), Cp_NBa = xi_NBaH2*Cp_NBaH2+
xi_NBaCO*Cp_NBaCO+xi_NBaC3H8*Cp_NBaC3H8+xi_NBaCH4*Cp_NBaCH4+xi_NBaH2O*
Cp_NBaH2O+xi_NBaCO2*Cp_NBaCO2+xi_NBaN2*Cp_NBaN2+xi_NBaO2*Cp_NBaO2, n_POXe = 
ni_POXeH2+ni_POXeCO+ni_POXeC3H8+ni_POXeCH4+ni_POXeH2O+ni_POXeCO2+ni_POXeN2+
ni_POXeO2, n_POXLuft = ni_POXLuftH2+ni_POXLuftCO+ni_POXLuftC3H8+ni_POXLuftCH4+
ni_POXLuftH2O+ni_POXLuftCO2+ni_POXLuftN2+ni_POXLuftO2, n_POXa = ni_POXaH2+
ni_POXaCO+ni_POXaC3H8+ni_POXaCH4+ni_POXaH2O+ni_POXaCO2+ni_POXaN2+ni_POXaO2, 
n_NBe = ni_NBeH2+ni_NBeCO+ni_NBeC3H8+ni_NBeCH4+ni_NBeH2O+ni_NBeCO2+ni_NBeN2+
ni_NBeO2, n_NBe2 = ni_NBe2H2+ni_NBe2CO+ni_NBe2C3H8+ni_NBe2CH4+ni_NBe2H2O+
ni_NBe2CO2+ni_NBe2N2+ni_NBe2O2, n_NBLuft = ni_NBLuftH2+ni_NBLuftCO+
ni_NBLuftC3H8+ni_NBLuftCH4+ni_NBLuftH2O+ni_NBLuftCO2+ni_NBLuftN2+ni_NBLuftO2,
n_NBa = ni_NBaH2+ni_NBaCO+ni_NBaC3H8+ni_NBaCH4+ni_NBaH2O+ni_NBaCO2+ni_NBaN2+
ni_NBaO2, N_POX = p_POXa*V_POX/R/T_POXa, Lambda_POXaH2 = .74525e-1+.388125e-4*
T_POXa-.30625e-7*T_POXa^2, Lambda_POXaCO = .1431e-1+.5215e-4*T_POXa, 
Lambda_POXaC3H8 = .1858e-2-.4698e-5*T_POXa+.2177e-6*T_POXa^2-.8409e-10*T_POXa^
3, Lambda_POXaCH4 = .6116e-1-.441e-4*T_POXa+.1525e-6*T_POXa^2, Lambda_POXaH2O
= -.1324e-1+.889e-4*T_POXa+.125e-7*T_POXa^2, Lambda_POXaCO2 = -.1557e-1+.10395\
e-3*T_POXa-.2125e-7*T_POXa^2, Lambda_POXaN2 = .5745e-2+.71925e-4*T_POXa-.15625\
e-7*T_POXa^2, Lambda_POXaO2 = .1294e-1+.5835e-4*T_POXa, Lambda_POX = (
Lambda_POXaH2*xi_POXaH2*Mi_H2^(1/3)+Lambda_POXaCO*xi_POXaCO*Mi_CO^(1/3)+
Lambda_POXaC3H8*xi_POXaC3H8*Mi_C3H8^(1/3)+Lambda_POXaCH4*xi_POXaCH4*Mi_CH4^(1/
3)+Lambda_POXaH2O*xi_POXaH2O*Mi_H2O^(1/3)+Lambda_POXaCO2*xi_POXaCO2*Mi_CO2^(1/
3)+Lambda_POXaN2*xi_POXaN2*Mi_N2^(1/3)+Lambda_POXaO2*xi_POXaO2*Mi_O2^(1/3))/(
xi_POXaH2*Mi_H2^(1/3)+xi_POXaCO*Mi_CO^(1/3)+xi_POXaC3H8*Mi_C3H8^(1/3)+
xi_POXaCH4*Mi_CH4^(1/3)+xi_POXaH2O*Mi_H2O^(1/3)+xi_POXaCO2*Mi_CO2^(1/3)+
xi_POXaN2*Mi_N2^(1/3)+xi_POXaO2*Mi_O2^(1/3)), Q_POXwand = A_POXwand*Nu_POXwand
*Lambda_POX/dh_POX*(T_POXwand-T_POXa), Q_NBwanda = A_NBwanda*Nu_NBwanda*
Lambda_POX/dh_POX*(T_NBwand-T_POXa), Q_POXNB = Q_POXwand+Q_NBwanda, N_NB = 
p_NBa*V_NB/R/T_NBa, Lambda_NBaH2 = .74525e-1+.388125e-4*T_NBa-.30625e-7*T_NBa^
2, Lambda_NBaCO = .1431e-1+.5215e-4*T_NBa, Lambda_NBaC3H8 = .1858e-2-.4698e-5*
T_NBa+.2177e-6*T_NBa^2-.8409e-10*T_NBa^3, Lambda_NBaCH4 = .6116e-1-.441e-4*
T_NBa+.1525e-6*T_NBa^2, Lambda_NBaH2O = -.1324e-1+.889e-4*T_NBa+.125e-7*T_NBa^
2, Lambda_NBaCO2 = -.1557e-1+.10395e-3*T_NBa-.2125e-7*T_NBa^2, Lambda_NBaN2 =
.5745e-2+.71925e-4*T_NBa-.15625e-7*T_NBa^2, Lambda_NBaO2 = .1294e-1+.5835e-4*
T_NBa, Lambda_NB = (Lambda_NBaH2*xi_NBaH2*Mi_H2^(1/3)+Lambda_NBaCO*xi_NBaCO*
Mi_CO^(1/3)+Lambda_NBaC3H8*xi_NBaC3H8*Mi_C3H8^(1/3)+Lambda_NBaCH4*xi_NBaCH4*
Mi_CH4^(1/3)+Lambda_NBaH2O*xi_NBaH2O*Mi_H2O^(1/3)+Lambda_NBaCO2*xi_NBaCO2*
Mi_CO2^(1/3)+Lambda_NBaN2*xi_NBaN2*Mi_N2^(1/3)+Lambda_NBaO2*xi_NBaO2*Mi_O2^(1/
3))/(xi_NBaH2*Mi_H2^(1/3)+xi_NBaCO*Mi_CO^(1/3)+xi_NBaC3H8*Mi_C3H8^(1/3)+
xi_NBaCH4*Mi_CH4^(1/3)+xi_NBaH2O*Mi_H2O^(1/3)+xi_NBaCO2*Mi_CO2^(1/3)+xi_NBaN2*
Mi_N2^(1/3)+xi_NBaO2*Mi_O2^(1/3)), Q_NBwandi = A_NBwandi*Nu_NBwandi*Lambda_NB/
dh_NB*(T_NBwand-T_NBa), Q_rad = Sigma*(T_NBwand^4-T_POXwand^4)/(1/A_NBwanda/
Epsilon_NBwanda+1/A_POXwand*(1/Epsilon_POXwand-1))],("Parameters")=[T_POXe = 
738.5015100518901, p_POXe = 101300, ni_POXeH2 = .267847489961e-3, ni_POXeCO =
.221582542707e-3, ni_POXeC3H8 = .37e-3, ni_POXeCH4 = .851480812e-6, ni_POXeH2O
= .33390675461e-3, ni_POXeCO2 = .230158881127e-3, ni_POXeN2 = .425652850783e-3
, ni_POXeO2 = 0, T_POXLuft = 298, p_POXLuft = 101300, ni_POXLuftH2 = 0, 
ni_POXLuftCO = 0, ni_POXLuftC3H8 = 0, ni_POXLuftCH4 = 0, ni_POXLuftH2O = 0, 
ni_POXLuftCO2 = 0, ni_POXLuftN2 = .1043928571428e-2, ni_POXLuftO2 = .2775e-3,
T_NBe = 1140.4504711677835, p_NBe = 101300, ni_NBeH2 = .328452691698e-3, 
ni_NBeCO = .271719487071e-3, ni_NBeC3H8 = 0, ni_NBeCH4 = .1044143309e-5, 
ni_NBeH2O = .409459025895e-3, ni_NBeCO2 = .282236372778e-3, ni_NBeN2 = .\
521964288667e-3, ni_NBeO2 = 0, T_NBLuft = 873, p_NBLuft = 101300, ni_NBLuftH2
= 0, ni_NBLuftCO = 0, ni_NBLuftC3H8 = 0, ni_NBLuftCH4 = 0, ni_NBLuftH2O = 0, 
ni_NBLuftCO2 = 0, ni_NBLuftN2 = .4175714285714e-2, ni_NBLuftO2 = .111e-2, 
T_NBe2 = 473, p_NBe2 = 101300, ni_NBe2H2 = 0, ni_NBe2CO = 0, ni_NBe2C3H8 = 0,
ni_NBe2CH4 = 0, ni_NBe2H2O = 0, ni_NBe2CO2 = 0, ni_NBe2N2 = 0, ni_NBe2O2 = 0,
xi_anaH2 = .180978033757249, xi_anaCO = .149717934261691, xi_anaC3H8 = 0, 
xi_anaCH4 = .575324872657e-3, xi_anaH2O = .225612672034038, xi_anaCO2 = .15551\
2757518103, xi_anaN2 = .287603277556262, xi_anaO2 = 0],("AlgVars")=[xi_POXaH2,
xi_POXaCO, xi_POXaCH4, xi_POXaH2O, xi_POXaCO2, xi_POXaN2],("DynVars")=[T_POXa,
T_NBa, T_POXwand, T_NBwand],("ODEs")=[`T_POXa'` = (Q_POXNB+ni_POXeH2*h_POXeH2+
ni_POXeCO*h_POXeCO+ni_POXeC3H8*h_POXeC3H8+ni_POXeCH4*h_POXeCH4+ni_POXeH2O*
h_POXeH2O+ni_POXeCO2*h_POXeCO2+ni_POXeN2*h_POXeN2+ni_POXeO2*h_POXeO2+
ni_POXLuftH2*h_POXLuftH2+ni_POXLuftCO*h_POXLuftCO+ni_POXLuftC3H8*h_POXLuftC3H8
+ni_POXLuftCH4*h_POXLuftCH4+ni_POXLuftH2O*h_POXLuftH2O+ni_POXLuftCO2*
h_POXLuftCO2+ni_POXLuftN2*h_POXLuftN2+ni_POXLuftO2*h_POXLuftO2-ni_POXaH2*
h_POXaH2-ni_POXaCO*h_POXaCO-ni_POXaC3H8*h_POXaC3H8-ni_POXaCH4*h_POXaCH4-
ni_POXaH2O*h_POXaH2O-ni_POXaCO2*h_POXaCO2-ni_POXaN2*h_POXaN2-ni_POXaO2*
h_POXaO2)/N_POX/Cp_POXa, `T_NBa'` = (ni_NBeH2*h_NBeH2+ni_NBeCO*h_NBeCO+
ni_NBeC3H8*h_NBeC3H8+ni_NBeCH4*h_NBeCH4+ni_NBeH2O*h_NBeH2O+ni_NBeCO2*h_NBeCO2+
ni_NBeN2*h_NBeN2+ni_NBeO2*h_NBeO2+ni_NBLuftH2*h_NBLuftH2+ni_NBLuftCO*
h_NBLuftCO+ni_NBLuftC3H8*h_NBLuftC3H8+ni_NBLuftCH4*h_NBLuftCH4+ni_NBLuftH2O*
h_NBLuftH2O+ni_NBLuftCO2*h_NBLuftCO2+ni_NBLuftN2*h_NBLuftN2+ni_NBLuftO2*
h_NBLuftO2+ni_NBe2H2*h_NBe2H2+ni_NBe2CO*h_NBe2CO+ni_NBe2C3H8*h_NBe2C3H8+
ni_NBe2CH4*h_NBe2CH4+ni_NBe2H2O*h_NBe2H2O+ni_NBe2CO2*h_NBe2CO2+ni_NBe2N2*
h_NBe2N2+ni_NBe2O2*h_NBe2O2-ni_NBaH2*h_NBaH2-ni_NBaCO*h_NBaCO-ni_NBaC3H8*
h_NBaC3H8-ni_NBaCH4*h_NBaCH4-ni_NBaH2O*h_NBaH2O-ni_NBaCO2*h_NBaCO2-ni_NBaN2*
h_NBaN2-ni_NBaO2*h_NBaO2+Q_NBwandi)/N_NB/Cp_NBa, `T_POXwand'` = (-Q_POXwand+
Q_rad)/Pi/(r4^2-r3^2)/LPOX/rho_POXwand/C_POXwand/por, `T_NBwand'` = -(
Q_NBwanda+Q_NBwandi+Q_rad)/Pi/(r2^2-r1^2)/LPOX/rho_NBwand/C_NBwand/por],("AEs"
)=[0 = 1-xi_POXaC3H8-xi_POXaH2-xi_POXaCO-xi_POXaCH4-xi_POXaH2O-xi_POXaCO2-
xi_POXaN2-xi_POXaO2, 0 = (8*xi_POXaC3H8+2*xi_POXaH2+4*xi_POXaCH4+2*xi_POXaH2O)
/(8*xi_POXmC3H8+2*xi_POXmH2+4*xi_POXmCH4+2*xi_POXmH2O)-(xi_POXaCO+xi_POXaH2O+2
*xi_POXaCO2+2*xi_POXaO2)/(xi_POXmCO+xi_POXmH2O+2*xi_POXmCO2+2*xi_POXmO2), 0 =
ni_POXLuftN2+ni_POXeN2-xi_POXaN2*n_POXa, 0 = (8*xi_POXaC3H8+2*xi_POXaH2+4*
xi_POXaCH4+2*xi_POXaH2O)/(8*xi_POXmC3H8+2*xi_POXmH2+4*xi_POXmCH4+2*xi_POXmH2O)
-(3*xi_POXaC3H8+xi_POXaCO+xi_POXaCH4+xi_POXaCO2)/(3*xi_POXmC3H8+xi_POXmCO+
xi_POXmCH4+xi_POXmCO2), 0 = -(G_H2+G_CO2-G_CO-G_H2O)/R/T_POXa-ln(xi_POXaH2*
xi_POXaCO2/xi_POXaCO/xi_POXaH2O), 0 = -(3*G_H2+G_CO-G_CH4-G_H2O)/R/T_POXa-2*ln
(p_POXa/p0)-ln(xi_POXaH2^3*p_POXa^2/p0^2*xi_POXaCO/xi_POXaCH4/xi_POXaH2O)]]);
getSys := proc () return eval(Sys) end proc; getSol := proc () return [
xi_POXaH2 = .368750587510893, xi_POXaCO = .284790826986732, xi_POXaCH4 = .\
22337314596529e-1, xi_POXaH2O = .12042268644267e-1, xi_POXaCO2 = .\
11972472427902e-1, xi_POXaN2 = .300106529833676, T_POXa = 1026.5889711019772,
T_NBa = 1166.1518337724226, T_POXwand = 1059.8353392024314, T_NBwand = 1076.26\
55452849394] end proc; end module; end module; Teymour1989b := module () 
export CSTR; CSTR := module () local Sys; export checkConsistency, getSys, 
getSolutions; Sys := table([("ExplicitAEs")=[yf = 2/3, yref = 5/9, Tc = 45, T
= 45*yT, negdeltaH = 21000, B = 1.4, alpha = 50/3, a3 = -.3495, a2 = -6.7530,
a1 = -.4407, f = .8, R = 1.987, MWs = 74.10, MWm = 86.05, phim = 1-phis, rhosf
= 777.2814025, rhomf = 918.5720, rhomg = 892.0200, rhof = rhomf*phim+rhosf*
phis, rhos = 74120/(91.878+.116*T), rhop = 1211-.8496*T, rhom = 958.4-1.3276*T
, rhofT = rhom*phim+rhos*phis, vp = 1-vm-vs, rho = rhom*vm+rhos*vs+rhop*vp, kd
= exp(ln(60)+34.99620124-30800/R/(T+273)), xt = (rhomf*phim-vm*rhomg)/(rhomf*
phim+MWm/MWs*phis*rhosf), lt = exp(a1*xt+a2*xt^2+a3*xt^3), ltrm = 417.3948000*
exp(-7569./R/(T+273.)), lp = 979.8*exp(-4869/R/(T+273))*15^(1/2), Cpf = (.470*
rhomf*phim+.716*rhosf*phis)/(rhomf*phim+rhosf*phis), Cp = (.470*rhom*vm+.716*
rhos*vs+rhop*vp*(.321425+.955e-3*T))/(rhom*vm+rhos*vs+rhop*vp), Q = 2^(1/2)*(f
*kd*zi*If/lt)^(1/2), Rm = (lp+ltrm)*Q*rhom/MWm*vm+2*f*kd*zi*If, qoutoverqin =
rhof/rhofT+theta*MWm*Rm*(1/rhop-1/rhom)],("Parameters")=[theta = 200.4617, If
= .316531e-1, phis = .6],("Ranges")=[vm = 0 .. 1, vs = 0 .. 1, zi = 0 .. 1, yT
= 0 .. 10, theta = 0 .. 300, If = 0 .. .1],("AlgVars")=[],("DynVars")=[vm, vs,
zi, yT],("ODEs")=[`vm'` = (1-phis)*rhomf/rhom-vm*qoutoverqin-theta*vm*(lp+ltrm
)*Q-2*theta*MWm/rhom*f*kd*zi*If, `vs'` = phis*rhosf/rhos-vs*qoutoverqin, `zi'`
= 1-theta*kd*zi-zi*qoutoverqin, `yT'` = rhof*Cpf/rho/Cp*(yf-yref)-rhof/rho*(yT
-yref)+negdeltaH/rho/Cp/Tc*theta*vm*lp*rhom/MWm*Q-theta*alpha*B/rho/Cp*(yT-1)]
,("AEs")=[]]); checkConsistency := proc () local Sol, Res, AbsRes, MaxAbsRes;
Sol := getSolutions()[1]; Res := evalf(subs(Sol,map(rhs,Sys["ODEs"]))); AbsRes
:= map(abs,Res); MaxAbsRes := max(op(AbsRes)); if not MaxAbsRes < .10e-5 then
error "error in ODEs larger than expected" end if; Res := evalf(subs(Sol,map(
lhs-rhs,Sys["ExplicitAEs"]))); AbsRes := map(abs,Res); MaxAbsRes := max(op(
AbsRes)); if not MaxAbsRes < .10e-5 then error 
"error in ExplicitAEs larger than expected" end if; return true end proc; 
getSys := proc () return eval(Sys) end proc; getSolutions := proc () local 
Sol1, Sol2, Sol3, Sol4, Sol5, Sol6, Equil; Sol1 := [vm = .1167175, vs = .64340\
86, zi = .9702816, yT = 1.250588, theta = 200.4617, If = .316531e-1, phis = .6
, yf = .6666666667, yref = .5555555556, Tc = 45., T = 56.276460, negdeltaH = 
21000., B = 1.4, alpha = 16.66666667, a3 = -.3495, a2 = -6.7530, a1 = -.4407,
f = .8, R = 1.987, MWs = 74.10, MWm = 86.05, phim = .4, rhosf = 777.2814025, 
rhomf = 918.5720, rhomg = 892.0200, rhof = 833.7976415, rhos = 753.2055744, 
rhop = 1163.187520, rhom = 883.6873717, rhofT = 805.3982933, vp = .2398739, 
rho = 866.7790518, kd = .3406505834e-3, xt = .2896722329, lt = .4951978566, 
ltrm = .3948005385e-2, lp = 2.224468994, Cpf = .6075954180, Cp = .5770130601,
Q = .5814096079e-2, Rm = .1554641634e-1, qoutoverqin = .9623413178]; Sol2 := [
vm = .1174607, vs = .643351, zi = .9642609, yT = 1.282066, theta = 175.2325, 
If = .316531e-1, phis = .6, yf = .6666666667, yref = .5555555556, Tc = 45., T
= 57.692970, negdeltaH = 21000., B = 1.4, alpha = 16.66666667, a3 = -.3495, a2
= -6.7530, a1 = -.4407, f = .8, R = 1.987, MWs = 74.10, MWm = 86.05, phim = .4
, rhosf = 777.2814025, rhomf = 918.5720, rhomg = 892.0200, rhof = 833.7976415,
rhos = 751.9499935, rhop = 1161.984053, rhom = 881.8068130, rhofT = 803.892721\
3, vp = .2391883, rho = 865.2784161, kd = .4167564465e-3, xt = .2889429224, lt
= .4968024699, ltrm = .4148571323e-2, lp = 2.296520430, Cpf = .6075954180, Cp
= .5775101293, Q = .6400511534e-2, Rm = .1774526010e-1, qoutoverqin = .9640345\
048, w1_Hopf[1] = .1922075917, w1_Hopf[2] = -.3111743177e-1, w1_Hopf[3] = .\
7017859892e-1, w1_Hopf[4] = .3145015097e-1, w2_Hopf[1] = -.1547885971, w2_Hopf
[2] = .2483523391e-1, w2_Hopf[3] = .2407872264e-2, w2_Hopf[4] = .9651899829, 
omega_Hopf = -3.832101806]; Sol3 := [vm = .1174607, vs = .643351, zi = .964260\
9, yT = 1.282066, theta = 175.2325, If = .316531e-1, phis = .6, yf = .66666666\
67, yref = .5555555556, Tc = 45., T = 57.692970, negdeltaH = 21000., B = 1.4,
alpha = 16.66666667, a3 = -.3495, a2 = -6.7530, a1 = -.4407, f = .8, R = 1.987
, MWs = 74.10, MWm = 86.05, phim = .4, rhosf = 777.2814025, rhomf = 918.5720,
rhomg = 892.0200, rhof = 833.7976415, rhos = 751.9499935, rhop = 1161.984053,
rhom = 881.8068130, rhofT = 803.8927213, vp = .2391883, rho = 865.2784161, kd
= .4167564465e-3, xt = .2889429224, lt = .4968024699, ltrm = .4148571323e-2, 
lp = 2.296520430, Cpf = .6075954180, Cp = .5775101293, Q = .6400511534e-2, Rm
= .1774526010e-1, qoutoverqin = .9640345048, w1_Hopf[1] = .2985380149e-1, 
w1_Hopf[2] = -.5446277708e-2, w1_Hopf[3] = .1719964531, w1_Hopf[4] = .27129224\
35, w2_Hopf[1] = -.8903305244, w2_Hopf[2] = .1436114129, w2_Hopf[3] = -.186190\
2261, w2_Hopf[4] = .2189003123, omega_Hopf = -3.832101806]; Sol4 := [vm = .\
4479726e-1, vs = .6561665, zi = .7439232, yT = 1.681804, theta = 81.01578, If
= .6117706e-1, phis = .6, yf = .6666666667, yref = .5555555556, Tc = 45., T =
75.681180, negdeltaH = 21000., B = 1.4, alpha = 16.66666667, a3 = -.3495, a2 =
-6.7530, a1 = -.4407, f = .8, R = 1.987, MWs = 74.10, MWm = 86.05, phim = .4,
rhosf = 777.2814025, rhomf = 918.5720, rhomg = 892.0200, rhof = 833.7976415, 
rhos = 736.3619774, rhop = 1146.701269, rhom = 857.9256654, rhofT = 784.987452\
6, vp = .29903624, rho = 864.5140164, kd = .4678228352e-2, xt = .3602483701, 
lt = .3494138181, ltrm = .7515910137e-2, lp = 3.365816406, Cpf = .6075954180,
Cp = .5772253550, Q = .3122403421e-1, Rm = .4738402245e-1, qoutoverqin = .9652\
150328, w1_Hopf[1] = .3380886713, w1_Hopf[2] = -.5797626875e-1, w1_Hopf[3] = .\
8262857783, w1_Hopf[4] = -.1056717986, w2_Hopf[1] = .2750254946, w2_Hopf[2] =
-.4570504704e-1, w2_Hopf[3] = -.1534885406, w2_Hopf[4] = -.2951836388, 
omega_Hopf = 5.651577654]; Sol5 := [If = .316531000000000035e-1, phis = .59999\
9999999999978, zi = .983006407814302152, yT = 1.14788577667081548, vs = .60926\
7654598932728, theta = 37.6197132983749754, vm = .332614224959867255, yf = .66\
66666667, yref = .5555555556, Tc = 45., T = 51.65485996, negdeltaH = 21000., B
= 1.4, alpha = 16.66666667, a3 = -.3495, a2 = -6.7530, a1 = -.4407, f = .8, R
= 1.987, MWs = 74.10, MWm = 86.05, phim = .4000000000, rhosf = 777.2814025, 
rhomf = 918.5720, rhomg = 892.0200, rhof = 833.7976415, rhos = 757.3314338, 
rhop = 1167.114031, rhom = 889.8230079, rhofT = 810.3280635, vp = .581181204e-\
1, rho = 825.2158104, kd = .1742902129e-3, xt = .7781035800e-1, lt = .92742682\
08, ltrm = .3348548156e-2, lp = 2.000866343, Cpf = .6075954180, Cp = .59939270\
69, Q = .3058742835e-2, Rm = .2109402742e-1, qoutoverqin = 1.010730619, w1_SN[
1] = .2468306963, w1_SN[2] = -.3697574414e-1, w1_SN[3] = .3172420182e-1, w1_SN
[4] = -.9678331348]; Sol6 := [If = .316531000000000035e-1, phis = .59999999999\
9999978, zi = .983006407814302152, yT = 1.14788577667081548, vs = .60926765459\
8932728, theta = 37.6197132983749754, vm = .332614224959867255, yf = .66666666\
67, yref = .5555555556, Tc = 45., T = 51.65485996, negdeltaH = 21000., B = 1.4
, alpha = 16.66666667, a3 = -.3495, a2 = -6.7530, a1 = -.4407, f = .8, R = 1.9\
87, MWs = 74.10, MWm = 86.05, phim = .4000000000, rhosf = 777.2814025, rhomf =
918.5720, rhomg = 892.0200, rhof = 833.7976415, rhos = 757.3314338, rhop = 
1167.114031, rhom = 889.8230079, rhofT = 810.3280635, vp = .581181204e-1, rho
= 825.2158104, kd = .1742902129e-3, xt = .7781035800e-1, lt = .9274268208, 
ltrm = .3348548156e-2, lp = 2.000866343, Cpf = .6075954180, Cp = .5993927069,
Q = .3058742835e-2, Rm = .2109402742e-1, qoutoverqin = 1.010730619, w1_SN[1] =
.9156583123, w1_SN[2] = -.1371674917, w1_SN[3] = .1176860477, w1_SN[4] = -.359\
0333245]; Equil := [Sol1, Sol2, Sol3, Sol4, Sol5, Sol6]; return Equil end proc
; end module; end module; end module;
