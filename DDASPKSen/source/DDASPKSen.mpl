DDASPKSen := module () local CreateAdiforFiles, createSharedObject, 
createExternalFunctions, EvalSen, prepareDirForDDASPK, RunAdifor, 
SubsIntoSimpleWrapper, Templates, InstanceCounter, UseAdifor; export 
appendSenParsToDAESys, init, CreateInstance, RestoreInstance, Tests; 
createSharedObject := proc (NameOfRes::string, resfun::procedure, NameOfJac::
string, jacfun::procedure, ddaspkDir::string, aDAEsys, InterfaceKeyword::name)
local oldDir, LinkerCommand, ddaspk, SharedObjectName, SharedObjectEnd, 
ObjectEnd, NumParsInDAEsys, NumSenParsInDAESys, MaxNumIndependentVars, 
NameOfResFor, NameOfResObj, NameOfJacFor, NameOfJacObj, status, soName, 
Filelist; SharedObjectEnd := ".so"; ObjectEnd := ".o"; NameOfResFor := cat(
NameOfRes,".f"); NameOfResObj := cat(NameOfRes,ObjectEnd); NameOfJacFor := cat
(NameOfJac,".f"); NameOfJacObj := cat(NameOfJac,ObjectEnd); oldDir := 
currentdir(); currentdir(ddaspkDir); if Aux:-FileOperations:-fileExists(
NameOfResFor) then currentdir(oldDir); error 
"routine %1 already exists in directory %2", NameOfResFor, ddaspkDir end if; 
if Aux:-FileOperations:-fileExists(NameOfJacFor) then currentdir(oldDir); 
error "routine %1 already exists in directory %2", NameOfJacFor, ddaspkDir end
if; NumParsInDAEsys := nops(aDAEsys["Parameters"]); NumSenParsInDAESys := nops
(aDAEsys["SenPars"]); MaxNumIndependentVars := 2*nops(aDAEsys["DynVars"])+2*
nops(aDAEsys["AlgVars"]); codegen[fortran](resfun,filename = NameOfResFor,mode
= double,precision = double); status := ssystem(cat(
"gfortran -shared -fPIC -O0 -c ",NameOfResFor," -o ",NameOfResObj)); if status
[1] <> 0 then system(cat("gfortran -fPIC -O0 -c ",NameOfResFor," -o ",
NameOfResObj)); error " in \"../maple/DDASPKSen/source/CreateSharedObject.txt\
\" command: gfortran -fPIC -O0 -c [NameOfResFor] -o [NameOfResObj]\"" end if;
codegen[fortran](jacfun,filename = "jac.f.tmp",mode = double,precision = 
double); status := ssystem(cat("sh ",_EnvModulesDir,
"DDASPKSen/ext_routines/linux64/shell_scripts/sedScript1 ",NameOfJacFor)); if
status[1] <> 0 then system(cat("sh ",_EnvModulesDir,
"DDASPKSen/ext_routines/linux64/shell_scripts/sedScript1 ",NameOfJacFor)); 
error 
" in \"../maple/DDASPKSen/source/CreateSharedObject.txt\" command: \"sedScript\
1\"" end if; status := ssystem("rm -f jac.f.tmp"); if status[1] <> 0 then 
system("rm -f jac.f.tmp"); error " in \"../maple/DDASPKSen/source/CreateShared\
Object.txt\" command: \"rm -f jac.f.tmp\"" end if; if UseAdifor = true then 
ADIFOR:-RunAD(NameOfRes,['y', 'yprime'],MaxNumIndependentVars,'delta','dense',
'[]',"j"); status := ssystem(cat("sh ",_EnvModulesDir,
"DDASPKSen/ext_routines/linux64/shell_scripts/sedScript2 ",NameOfResFor)); if
status[1] <> 0 then system(cat("sh ",_EnvModulesDir,
"DDASPKSen/ext_routines/linux64/shell_scripts/sedScript2 ",NameOfResFor)); 
error 
" in \"../maple/DDASPKSen/source/createSharedObject.txt\" command: \"sedScript\
2\"" end if; status := ssystem(cat("rm -f j_",NameOfResFor)); if status[1] <>
0 then error " in \"../maple/DDASPKSen/source/createSharedObject.txt\" command\
: \"rm -f j_[NameOfResFor]\""; system(cat("rm -f j_",NameOfResFor)) end if end
if; status := ssystem(cat("gfortran -fPIC -O0 -c ",NameOfJacFor," -o ",
NameOfJacObj)); if status[1] <> 0 then system(cat("gfortran -fPIC -O0 -c ",
NameOfJacFor," -o ",NameOfJacObj)); error " in \"../maple/DDASPKSen/source/cre\
ateSharedObject.txt\" command: \"gfortran -fPIC -O0 -c [NameOfJacFor] -o [Name\
OfJacObj]\"" end if; if UseAdifor = true then status := ssystem(cat(
"gfortran -fPIC -O0 -c j",NameOfResFor," -o j",NameOfResObj)); if status[1] <>
0 then system(cat("gfortran -fPIC -O0 -c j",NameOfResFor," -o j",NameOfResObj)
); error " in \"../maple/DDASPKSen/source/createSharedObject.txt\" command: \"\
gfortran -fPIC -O0 -c j[NameOfResFor] -o j[NameOfResObj]\"" end if end if; if
EvalSen = true then ADIFOR:-RunAD(NameOfRes,['y', 'yprime', 'senpar'],
MaxNumIndependentVars+NumSenParsInDAESys,'delta'); status := ssystem(cat("sh "
,_EnvModulesDir,"DDASPKSen/ext_routines/linux64/shell_scripts/sedScript3 ",
NameOfResFor)); if status[1] <> 0 then system(cat("sh ",_EnvModulesDir,
"DDASPKSen/ext_routines/linux64/shell_scripts/sedScript3 ",NameOfResFor)); 
error 
" in \"../maple/DDASPKSen/source/createSharedObject.txt\" command: \"sedScript\
3\"" end if; status := ssystem(cat("rm -f g_",NameOfResFor)); if status[1] <>
0 then system(cat("rm -f g_",NameOfResFor)); error " in \"../maple/DDASPKSen/s\
ource/createSharedObject.txt\" command: \"rm -f g_[NameOfResFor]\"" end if; 
status := ssystem(cat("gfortran -fPIC -O0 -c g",NameOfResFor," -o g",
NameOfResObj)); if status[1] <> 0 then system(cat("gfortran -fPIC -O0 -c g",
NameOfResFor," -o g",NameOfResObj)); error " in \"../maple/DDASPKSen/source/cr\
eateSharedObject.txt\" command: \"gfortran -fPIC -O0 -c g[NameOfResFor] -o g[N\
ameOfResObj]\"" end if end if; if InterfaceKeyword = ('BuildInterface') then 
SharedObjectName := cat("cwrap_ddaspk",InstanceCounter,SharedObjectEnd); if 
UseAdifor = false and EvalSen = false then LinkerCommand := cat(
"gcc -shared -Xlinker -Bsymbolic -m64 ",_EnvModulesDir,
"DDASPKSen/ext_routines/f2c/shared_obj/ddaspk.so ",_EnvModulesDir,
"DDASPKSen/ext_routines/f2c/shared_obj/adf_dummy.o ",_EnvModulesDir,
"DDASPKSen/ext_routines/f2c/cwrap/g_res_dummy.o ",_EnvModulesDir,
"DDASPKSen/ext_routines/f2c/cwrap/psol_dummy.o ",_EnvModulesDir,
"DDASPKSen/ext_routines/f2c/cwrap/cwrap_ddaspk.o ","./",NameOfResObj," ./",
NameOfJacObj," ","-lm -o ",SharedObjectName) elif UseAdifor = true and EvalSen
= false then LinkerCommand := cat("sh ",_EnvModulesDir,"DDASPKSen/ext_routines\
/linux64/shell_scripts/linkerCommand_UseAdifor_True_EvalSen_False ",
_EnvModulesDir," ",NameOfResObj," ",NameOfJacObj," ",SharedObjectName) elif 
UseAdifor = true and EvalSen = true then LinkerCommand := cat("sh ",
_EnvModulesDir,"DDASPKSen/ext_routines/linux64/shell_scripts/linkerCommand_Use\
Adifor_True_EvalSen_True ",_EnvModulesDir," ",NameOfResObj," ",NameOfJacObj,
" ",SharedObjectName) end if; status := ssystem(eval(LinkerCommand)); if 
status[1] <> 0 then printf(
"------------------------------------------------------------------\n"); 
printf("%s\n",LinkerCommand); printf(
"------------------------------------------------------------------\n"); error
" in \"../maple/DDASPKSen/source/CreateSharedObject.txt\" linking failed!" end
if; print(currentdir()); soName := cat(_EnvWorkingDir,"/",ddaspkDir,"/",
SharedObjectName); ddaspk := define_external('cwrapper_ddaspk',NEQ::integer[4]
,T::REF(float[8]),Y::ARRAY(1 .. NEQ,float[8]),YPRIME::ARRAY(1 .. NEQ,float[8])
,tout::float[8],INFO::ARRAY(1 .. 30,integer[4]),RTOL::REF(float[8]),ATOL::REF(
float[8]),IDID::REF(integer[4]),RWORK::ARRAY(1 .. LENW,float[8]),LENW::integer
[4],IWORK::ARRAY(1 .. LENIW,integer[4]),LENIW::integer[4],RPAR::ARRAY(1 .. 
NumParsInDAEsys,float[8]),IPAR::integer[4],SENPAR::ARRAY(1 .. 
NumSenParsInDAESys,float[8]),LIB = soName); Filelist := [convert(NameOfRes,
name) = convert(NameOfResObj,name), convert(NameOfJac,name) = convert(
NameOfJacObj,name)]; if UseAdifor = true then if EvalSen = true then Filelist
:= [op(Filelist), cat(g,NameOfRes) = cat(g,NameOfResObj)] end if; Filelist :=
[op(Filelist), cat(j,NameOfRes) = cat(j,NameOfResObj)] end if; save Filelist,
"Filelist.m"; currentdir(oldDir); return ddaspk else Filelist := [convert(
NameOfRes,name) = convert(NameOfResObj,name), convert(NameOfJac,name) = 
convert(NameOfJacObj,name)]; if UseAdifor = true then if EvalSen = true then 
Filelist := [op(Filelist), cat(g,NameOfRes) = cat(g,NameOfResObj)] end if; 
Filelist := [op(Filelist), cat(j,NameOfRes) = cat(j,NameOfResObj)] end if; 
save Filelist, "Filelist.m"; currentdir(oldDir) end if end proc; 
createExternalFunctions := proc (aDAESys::DAESys, NameOfRes::name, NameOfJac::
name) local i1, ListListOfJacElements, Eqs, EqsAfterSubs, 
NameOfDerivativeRoutine, RhsODEs, RhsAEs, NumExplicitAEs, 
ExplicitAEsPlusDummyForResTemplate, NumODEs, NumAEs, NumDynVars, NumAlgVars, 
NumPars, NumSenPars, UnitNumber; UnitNumber := 6; NameOfDerivativeRoutine := 
cat(j,NameOfRes); NumODEs := nops(aDAESys["ODEs"]); NumAEs := nops(aDAESys[
"AEs"]); NumExplicitAEs := nops(aDAESys["ExplicitAEs"]); NumDynVars := nops(
aDAESys["DynVars"]); NumAlgVars := nops(aDAESys["AlgVars"]); NumPars := max(1,
nops(aDAESys["Parameters"])); NumSenPars := max(1,nops(aDAESys["SenPars"])); 
RhsODEs := map(rhs,aDAESys["ODEs"]); RhsODEs := convert(RhsODEs,array); for i1
to NumODEs do RhsODEs[i1] := RhsODEs[i1]-yprime[i1] end do; RhsODEs := convert
(RhsODEs,list); RhsODEs := map(evalf,RhsODEs); RhsAEs := map(rhs,aDAESys["AEs"
]); RhsAEs := map(evalf,RhsAEs); Eqs := [op(RhsODEs), op(RhsAEs)]; if aDAESys[
"ExplicitAEs"] = [] then ExplicitAEsPlusDummyForResTemplate := [0.] else 
ExplicitAEsPlusDummyForResTemplate := [op(map(rhs,aDAESys["ExplicitAEs"])), 0.
] end if; ExplicitAEsPlusDummyForResTemplate := map(evalf,
ExplicitAEsPlusDummyForResTemplate); NameOfRes := subs(NumODEsToBeSubs = 
NumODEs,NumAEsToBeSubs = NumAEs,NumExplicitAEsToBeSubs = NumExplicitAEs,
NumDynVarsToBeSubs = NumDynVars,NumAlgVarsToBeSubs = NumAlgVars,
NumParsToBeSubs = NumPars,NumSenParsToBeSubs = NumSenPars,ExplicitAEsToBeSubs
= ExplicitAEsPlusDummyForResTemplate,ResToBeSubs = Eqs,eval(Templates:-
resTemplate)); if UseAdifor = true then NameOfJac := subs(NumODEsToBeSubs = 
NumODEs,NumAEsToBeSubs = NumAEs,NumDynVarsToBeSubs = NumDynVars,
NumAlgVarsToBeSubs = NumAlgVars,NumParsToBeSubs = NumPars,NumSenParsToBeSubs =
NumSenPars,UnitNumberToBeSubs = UnitNumber,NameOfDerivateRoutineToBeSubs = 
NameOfDerivativeRoutine,eval(Templates:-jacAdiforTemplate)) else EqsAfterSubs
:= [seq(Aux:-ListOperations:-subsEqnListIntoEqn(aDAESys["ExplicitAEs"],0 = Eqs
[i1]),i1 = 1 .. NumODEs+NumAEs)]; EqsAfterSubs := map(rhs,EqsAfterSubs); 
ListListOfJacElements := [seq([seq(diff(EqsAfterSubs[i1],y[i2])+cj*diff(
EqsAfterSubs[i1],yprime[i2]),i2 = 1 .. NumDynVars+NumAlgVars)],i1 = 1 .. 
NumODEs+NumAEs)]; NameOfJac := subs(NumODEsToBeSubs = NumODEs,NumAEsToBeSubs =
NumAEs,NumDynVarsToBeSubs = NumDynVars,NumAlgVarsToBeSubs = NumAlgVars,
NumParsToBeSubs = NumPars,NumSenParsToBeSubs = NumSenPars,JacToBeSubs = 
ListListOfJacElements,eval(Templates:-jacTemplate)) end if; RETURN() end proc;
prepareDirForDDASPK := proc (NameOfDir::string) if not Aux:-FileOperations:-
dirExists(NameOfDir) then mkdir(NameOfDir) end if; if Aux:-FileOperations:-
fileExists(cat(NameOfDir,"/mwrap*")) then return false else return true end if
end proc; SubsIntoSimpleWrapper := proc (aDAESys::DAESys, DDASPKtoBeSubs::
procedure) local NY, base, LRWtoBeSubs, LIWtoBeSubs, INFO5toBeSubs, 
INFO11toBeSubs, INFO19toBeSubs, INFO20toBeSubs, INFO22toBeSubs, INFO23toBeSubs
, INFO24toBeSubs, INFO25toBeSubs, NEQtoBeSubs, SimpleWrapper, 
ParsInDAEsysToBeSubs, SenParsInDAESysToBeSubs, IWORK38, 
ActInputPositionToBeSubs, SysInputToBeSubs, NameInputParameterToBeSubs; 
ActInputPositionToBeSubs := 0; SysInputToBeSubs := []; 
NameInputParameterToBeSubs := 0; NEQtoBeSubs := (nops(aDAESys["ODEs"])+nops(
aDAESys["AEs"]))*(nops(aDAESys["SenPars"])+1); NY := NEQtoBeSubs/(nops(aDAESys
["SenPars"])+1); ParsInDAEsysToBeSubs := map(rhs,aDAESys["Parameters"]); 
ParsInDAEsysToBeSubs := map(evalf,ParsInDAEsysToBeSubs); 
SenParsInDAESysToBeSubs := map(rhs,aDAESys["SenPars"]); 
SenParsInDAESysToBeSubs := map(evalf,SenParsInDAESysToBeSubs); if UseAdifor =
true then INFO5toBeSubs := 1 else INFO5toBeSubs := 1 end if; if EvalSen = true
then INFO19toBeSubs := nops(aDAESys["SenPars"]); INFO20toBeSubs := 3; 
INFO22toBeSubs := nops(aDAESys["SenPars"]); INFO23toBeSubs := 0; 
INFO24toBeSubs := 0; INFO25toBeSubs := 0 else INFO19toBeSubs := 0; 
INFO20toBeSubs := 0; INFO22toBeSubs := 0; INFO23toBeSubs := 0; INFO24toBeSubs
:= 0; INFO25toBeSubs := 0 end if; base := 50+9*NEQtoBeSubs; base := base+
NEQtoBeSubs; base := base+NY*NY; base := base+INFO19toBeSubs*(2*NY+max(NY,
INFO24toBeSubs)+INFO22toBeSubs); LRWtoBeSubs := base; base := 40+NY; 
LIWtoBeSubs := base; if not aDAESys["AEs"] = [] then INFO11toBeSubs := 1; 
LIWtoBeSubs := 40+2*NEQtoBeSubs end if; SimpleWrapper := module () local 
setActInputPosition, isApostrophe; export init, doIntStep, Sys, getY, 
getYPRIME, getT, getSys, getParameters, getVariables, getMaxResidueODEs, 
getSenParameters, getSenPars, getSensitivities, getS, getSPRIME, 
getMaxResidueExplicitAEs, getExplicitAEs, evalExplicitAEsInDAESys, 
setParameters, setVars, setVariables, setSenPars, setSenParameters, 
setInitialTime, setEndTime, runInt, runIntWithDynConstraint, setMaxNumSteps, 
checkInputForCalculateIntResults, checkInputForIntegrateAtSpecificTimePoints,
createCFunction, calculateIntResults, setInput, getInputCurve, 
getNameInputParameter, getActInputPosition, integrate, 
integrateAtSpecificTimePoints, setSpecificParameters, NEQ, LIW, LRW, DDASPK, N
, setPars, getPars, rpar, ipar, senpar, idid, t_vec, y_vec, yprime_vec, 
tout_vec, INFO_vec, RWORK_vec, IWORK_vec, rtol_vec, atol_vec, rpar_vec, 
ipar_vec, calculateYPRIME, senpar_vec, INFO5, INFO11, INFO19, INFO20, INFO22,
INFO23, INFO24, INFO25, NumParsInDAEsys, ParsInDAEsys, NumSenParsInDAESys, 
SenParsInDAESys, MaxNumSteps, VarNames, NumVars, data, item, i1, DdaspkStatus,
CurrentPoint, i2, ActInputPosition, SysInput, NameInputParameter; 
ActInputPosition := ActInputPositionToBeSubs; SysInput := SysInputToBeSubs; 
NameInputParameter := NameInputParameterToBeSubs; INFO5 := INFO5toBeSubs; 
INFO19 := INFO19toBeSubs; INFO20 := INFO20toBeSubs; INFO22 := INFO22toBeSubs;
INFO23 := INFO23toBeSubs; INFO24 := INFO24toBeSubs; INFO25 := INFO25toBeSubs;
NEQ := NEQtoBeSubs; LRW := LRWtoBeSubs; LIW := LIWtoBeSubs; DDASPK := 
DDASPKtoBeSubs; Sys := aDAESys; ParsInDAEsys := ParsInDAEsysToBeSubs; 
SenParsInDAESys := SenParsInDAESysToBeSubs; NumParsInDAEsys := nops(
ParsInDAEsys); NumSenParsInDAESys := nops(SenParsInDAESys); idid := 0; N := 30
; INFO_vec := Vector(1 .. N,fill = 0,datatype = integer[4]); INFO_vec[3] := 1;
INFO_vec[5] := INFO5; INFO_vec[19] := INFO19; INFO_vec[20] := INFO20; INFO_vec
[22] := INFO22; INFO_vec[23] := INFO23; INFO_vec[24] := INFO24; INFO_vec[25] 
:= INFO25; RWORK_vec := Vector(1 .. LRW,fill = 0,datatype = float[8]); 
IWORK_vec := Vector(1 .. LIW,fill = 0,datatype = integer[4]); rtol_vec := .10e\
-5; atol_vec := .10e-5; rpar_vec := Vector(1 .. NumParsInDAEsys,map(evalf,
ParsInDAEsys),datatype = float[8]); ipar_vec := 0; senpar_vec := Vector(1 .. 
NumSenParsInDAESys,map(evalf,SenParsInDAESys),datatype = float[8]); 
MaxNumSteps := 100; t_vec := 0; tout_vec := 1; if not aDAESys["AEs"] = [] then
INFO11 := INFO11toBeSubs; INFO_vec[11] := INFO11; for i1 to nops(Sys["DynVars"
]) do IWORK_vec[40+i1] := 1 end do; for i1 to nops(Sys["AlgVars"]) do 
IWORK_vec[40+nops(Sys["DynVars"])+i1] := -1 end do end if; calculateYPRIME :=
proc () local SubsList, yprime, aSys, dgdp, dgdy, dSdt, DimOldIWORK, i1, 
NewIWORK_vec; SubsList := [seq(Sys["DynVars"][i1] = y_vec[i1],i1 = 1 .. nops(
Sys["DynVars"])), seq(Sys["AlgVars"][i1] = y_vec[nops(Sys["DynVars"])+i1],i1 =
1 .. nops(Sys["AlgVars"])), seq(lhs(Sys["Parameters"][i1]) = rpar_vec[i1],i1 =
1 .. nops(Sys["Parameters"])), seq(lhs(Sys["SenPars"][i1]) = senpar_vec[i1],i1
= 1 .. nops(Sys["SenPars"]))]; yprime_vec := Vector(1 .. NEQ,i1 -> 0,datatype
= float[8]); if Sys["AEs"] = [] then try yprime := Aux:-SystemClasses:-
evalODEsInDAESys(Sys,SubsList); for i1 to nops(yprime) do yprime_vec[i1] := 
yprime[i1] end do; aSys := Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(Sys);
dgdp := [seq([seq(diff(rhs(aSys["ODEs"][i1]),lhs(Sys["SenPars"][i2])),i2 = 1 
.. nops(Sys["SenPars"]))],i1 = 1 .. nops(aSys["ODEs"]))]; dgdy := [seq(diff(
rhs(aSys["ODEs"][i1]),aSys["DynVars"][i1]),i1 = 1 .. nops(aSys["ODEs"]))]; 
dSdt := [seq([seq(dgdy[i2]*y_vec[nops(Sys["DynVars"])+(i1-1)*nops(Sys[
"DynVars"])+i2]+dgdp[i2,i1],i2 = 1 .. nops(Sys["DynVars"]))],i1 = 1 .. nops(
Sys["SenPars"]))]; dSdt := evalf(subs([op(aSys["Parameters"]), op(Sys[
"SenPars"]), seq(Sys["DynVars"][i1] = y_vec[i1],i1 = 1 .. nops(aSys["DynVars"]
))],dSdt)); for i1 to nops(Sys["SenPars"]) do for i2 to nops(Sys["DynVars"]) 
do yprime_vec[nops(Sys["DynVars"])+(i1-1)*nops(Sys["DynVars"])+i2] := dSdt[i1]
[i2] end do end do catch: printf("\n WARNING: Initialization failed \n"); 
INFO_vec[11] := 1; LIW := 40+2*nops(Sys["DynVars"]); NewIWORK_vec := Vector(1
.. LIW,fill = 0,datatype = integer[4]); DimOldIWORK := LinearAlgebra[Dimension
](IWORK_vec); for i1 to DimOldIWORK do NewIWORK_vec[i1] := IWORK_vec[i1] end 
do; IWORK_vec := NewIWORK_vec; for i1 to nops(Sys["DynVars"]) do IWORK_vec[40+
i1] := 1 end do end try end if; RETURN() end proc; 
checkInputForCalculateIntResults := proc (timePoints::list(numeric), initTime
::numeric) local timeIsOrdered; if not initTime < timePoints[1] then error "th\
e first TOUT value in timePoints (%1) has to be smaller than initial time (%2)\
", timePoints[1], initTime end if; timeIsOrdered := ListTools[Sorted](
timePoints,verify,less_than); if not timeIsOrdered then error 
"Input timePoints has to be a strictly ascending list. You have inserted: %1",
timePoints end if end proc; checkInputForIntegrateAtSpecificTimePoints := proc
(timePoints::list(numeric), initTime::numeric) local timeIsOrdered; if not 
initTime < timePoints[1] then error "the first TOUT value in timePoints (%1) h\
as to be smaller than initial time (%2)", timePoints[1], initTime end if; 
timeIsOrdered := ListTools[Sorted](timePoints,verify,less_than); if not 
timeIsOrdered then error 
"Input timePoints has to be a strictly ascending list. You have inserted: %1",
timePoints end if end proc; doIntStep := proc () DDASPK(NEQ,'t_vec',y_vec,
yprime_vec,tout_vec,INFO_vec,'rtol_vec','atol_vec','idid',RWORK_vec,LRW,
IWORK_vec,LIW,rpar_vec,ipar_vec,senpar_vec); RETURN(idid) end proc; 
evalExplicitAEsInDAESys := proc () local SubsList; SubsList := Aux:-
SystemClasses:-evalExplicitAEsInDAESys(getSys(),[op(getParameters()), op(
getSenParameters()), op(getVariables())]); return SubsList end proc; 
getExplicitAEs := proc () local ParsAndVars, SubsList, EAEs, item, NewItem; 
ParsAndVars := [op(getParameters()), op(getSenParameters()), op(getVariables()
)]; SubsList := []; EAEs := Sys["ExplicitAEs"]; for item in EAEs while true do
NewItem := subs(ParsAndVars,item); NewItem := Aux:-ListOperations:-
subsEqnListIntoEqn(SubsList,NewItem); SubsList := [op(SubsList), NewItem] end
do; return SubsList end proc; getParameters := proc () local Pars, NamesOfPars
, SubsList; Pars := getPars(); NamesOfPars := Aux:-ListOperations:-
getListOfLhsIn(Sys["Parameters"]); SubsList := [seq(NamesOfPars[i1] = Pars[i1]
,i1 = 1 .. nops(Sys["Parameters"]))]; return SubsList end proc; getPars := 
proc () return convert(rpar_vec,list) end proc; getSenParameters := proc () 
local SenParameters, NamesOfSenPars, SubsList, NamesOfPars; SenParameters := 
getSenPars(); NamesOfPars := Aux:-ListOperations:-getListOfLhsIn(Sys["SenPars"
]); SubsList := [seq(NamesOfPars[i1] = SenParameters[i1],i1 = 1 .. nops(Sys[
"SenPars"]))]; return SubsList end proc; getSenPars := proc () return convert(
senpar_vec,list) end proc; getMaxResidueODEs := proc () local Res, AbsRes, 
MaxAbsRes, SolIncludingEAEs, PosOfMax, i1; SolIncludingEAEs := Aux:-
SystemClasses:-evalExplicitAEsInDAESys(getSys(),[op(getParameters()), op(
getSenParameters()), op(getVariables())]); Res := evalf(subs(SolIncludingEAEs,
map(rhs,getSys()["ODEs"]))); AbsRes := map(abs,Res); MaxAbsRes := 0; PosOfMax
:= 1; for i1 to nops(AbsRes) do if MaxAbsRes < AbsRes[i1] then MaxAbsRes := 
AbsRes[i1]; PosOfMax := i1 end if end do; return Sys["DynVars"][PosOfMax] = 
MaxAbsRes end proc; getS := proc () local Y, NumDynVars, NumAlgVars, 
NumSenPars; Y := getY(); NumDynVars := nops(Sys["DynVars"]); NumAlgVars := 
nops(Sys["AlgVars"]); NumSenPars := nops(Sys["SenPars"]); return Y[NumDynVars+
NumAlgVars+1 .. -1] end proc; getSensitivities := proc () local Y, SubsTable,
i1, i2, Sens, AlgSensitivities, DynSensitivities, NumDynVars, NumAlgVars, 
NumSenPars, SenDynNames, SenAlgNames; Y := getY(); SubsTable := table(); 
NumDynVars := nops(Sys["DynVars"]); NumAlgVars := nops(Sys["AlgVars"]); 
NumSenPars := nops(Sys["SenPars"]); Sens := getS(); SenDynNames := [seq(seq(
convert(cat(d,convert(Sys["DynVars"][i2],string),d,convert(lhs(Sys["SenPars"][
i1]),string)),name),i2 = 1 .. nops(Sys["DynVars"])),i1 = 1 .. nops(Sys[
"SenPars"]))]; SenAlgNames := [seq(seq(convert(cat(d,convert(Sys["AlgVars"][i2
],string),d,convert(lhs(Sys["SenPars"][i1]),string)),name),i2 = 1 .. nops(Sys[
"AlgVars"])),i1 = 1 .. nops(Sys["SenPars"]))]; DynSensitivities := seq(seq(
SenDynNames[i1+(i2-1)*nops(Sys["DynVars"])] = Sens[i1+(i2-1)*(nops(Sys[
"DynVars"])+nops(Sys["AlgVars"]))],i1 = 1 .. nops(Sys["DynVars"])),i2 = 1 .. 
nops(Sys["SenPars"])); AlgSensitivities := seq(seq(SenAlgNames[i1+(i2-1)*nops(
Sys["AlgVars"])] = Sens[i1+nops(Sys["DynVars"])+(i2-1)*(nops(Sys["DynVars"])+
nops(Sys["AlgVars"]))],i1 = 1 .. nops(Sys["AlgVars"])),i2 = 1 .. nops(Sys[
"SenPars"])); return [DynSensitivities, AlgSensitivities] end proc; getSPRIME
:= proc () local YPrime, NumDynVars, NumAlgVars, NumSenPars; YPrime := 
getYPRIME(); NumDynVars := nops(Sys["DynVars"]); NumAlgVars := nops(Sys[
"AlgVars"]); NumSenPars := nops(Sys["SenPars"]); return YPrime[NumDynVars+
NumAlgVars+1 .. -1] end proc; getMaxResidueExplicitAEs := proc () local Res, 
AbsRes, MaxAbsRes, SolIncludingEAEs; SolIncludingEAEs := Aux:-SystemClasses:-
evalExplicitAEsInDAESys(getSys(),[op(getParameters()), op(getSenParameters()),
op(getVariables())]); Res := evalf(subs(SolIncludingEAEs,map(lhs-rhs,getSys()[
"ExplicitAEs"]))); AbsRes := map(abs,Res); MaxAbsRes := max(op(AbsRes)); 
return MaxAbsRes end proc; getSys := proc () return eval(Sys) end proc; getT 
:= proc () RETURN(t_vec) end proc; getVariables := proc () local Y, SubsList;
Y := getY(); SubsList := [seq(Sys["DynVars"][i1] = Y[i1],i1 = 1 .. nops(Sys[
"DynVars"])), seq(Sys["AlgVars"][i1-nops(Sys["DynVars"])] = Y[i1],i1 = nops(
Sys["DynVars"])+1 .. nops(Sys["DynVars"])+nops(Sys["AlgVars"]))]; return 
SubsList end proc; getY := proc () RETURN(convert(y_vec,list)) end proc; 
getYPRIME := proc () RETURN(convert(yprime_vec,list)) end proc; init := proc (
) calculateYPRIME(); INFO_vec[1] := 0; INFO_vec[3] := 1; idid := doIntStep();
RETURN(idid) end proc; integrate := proc () local NumVarsInDAESys, Pars, 
ValuesExplicitVars, OtherData, SenNames, Curve, printInfo; printInfo := true;
if nargs = 1 then printInfo := args[1]; if not hastype(printInfo,boolean) then
error 
"Optional input argument has to be of type boolean. You have inserted %1.", 
printInfo end if elif 1 < nargs then error 
"Too many optional arguments. Only one optional argument is allowed!" end if;
SenNames := map(op,[seq([seq(convert(cat(d,convert(Sys["DynVars"][i2],string),
d,convert(lhs(Sys["SenPars"][i1]),string)),name),i2 = 1 .. nops(Sys["DynVars"]
)), seq(convert(cat(d,convert(Sys["AlgVars"][i2],string),d,convert(lhs(Sys[
"SenPars"][i1]),string)),name),i2 = 1 .. nops(Sys["AlgVars"]))],i1 = 1 .. nops
(Sys["SenPars"]))]); VarNames := [op(Sys["DynVars"]), op(Sys["AlgVars"]), op(
SenNames), op(map(lhs,Sys["ExplicitAEs"]))]; NumVars := nops(VarNames); data 
:= table(); if member(t,VarNames) then error 
"name t must not be a variable names" end if; for i1 to MaxNumSteps do if 
printInfo then if `mod`(i1,100) = 0 then printf(
"done with %d out of %d steps\n",i1,MaxNumSteps) end if end if; DdaspkStatus 
:= doIntStep(); ValuesExplicitVars := map(evalf,map(rhs,getExplicitAEs())); 
ValuesExplicitVars := subs(t = getT(),ValuesExplicitVars); CurrentPoint := [op
(getY()), op(ValuesExplicitVars)]; if not DdaspkStatus = 1 then for i2 to nops
(VarNames) do item := VarNames[i2]; data[item] := CurrentPoint[i2] end do; 
data[t] := getT(); break end if end do; if printInfo then print(
"Integrate finishes with DDASPK status ",DdaspkStatus) end if; return [data, 
DdaspkStatus] end proc; integrateAtSpecificTimePoints := proc (timePoints::
list(numeric)) local intiY, initTime, initPoint, timeIsOrdered, senNames, 
varNames, resultTable, results, nrTimePoints, actData, actName, actTOUT, 
actResult, statusInt, lastPosition, otherData, curve, valuesExplicitVars, i, 
inputTimes, inputValues, intAndInputTimes, k, initY, 
positionJustSmallerThanEnd, entriesSmallerLastTimePoint, printInfo, endTime, 
sysVarNames, timeTable, inputTable, integrationTime, actInputIndex, evalbList,
actInputValue, allTimePoints, integrationCounter, numberOfVariables; printInfo
:= true; if nargs = 2 then printInfo := args[2]; if not hastype(printInfo,
boolean) then error 
"Optional input argument has to be of type boolean. You have inserted %1.", 
printInfo end if elif 2 < nargs then error 
"Too many optional arguments. Only one optional argument is allowed!" end if;
nrTimePoints := nops(timePoints); endTime := timePoints[-1]; initY := getY();
initTime := getT(); valuesExplicitVars := map(evalf,map(rhs,getExplicitAEs()))
; checkInputForCalculateIntResults(timePoints,initTime); if not type(initY,
list(numeric)) then error 
"Initial values for variables have not been set. Call SetVariables()." end if;
calculateYPRIME(); INFO_vec[1] := 0; senNames := map(op,[seq([seq(convert(cat(
d,convert(Sys["DynVars"][i2],string),d,convert(lhs(Sys["SenPars"][i1]),string)
),name),i2 = 1 .. nops(Sys["DynVars"])), seq(convert(cat(d,convert(Sys[
"AlgVars"][i2],string),d,convert(lhs(Sys["SenPars"][i1]),string)),name),i2 = 1
.. nops(Sys["AlgVars"]))],i1 = 1 .. nops(Sys["SenPars"]))]); sysVarNames := [
op(Sys["DynVars"]), op(Sys["AlgVars"])]; varNames := [op(sysVarNames), op(
senNames), op(map(lhs,Sys["ExplicitAEs"]))]; resultTable := table(); 
resultTable[time] := array(1 .. nrTimePoints+1); resultTable[time][1] := 
initTime; initPoint := [op(initY), op(valuesExplicitVars)]; for i to nops(
varNames) do actName := varNames[i]; resultTable[actName] := array(1 .. 
nrTimePoints+1); resultTable[actName][1] := initPoint[i] end do; timeTable :=
table([seq(evalf(timePoints[i]) = 1,i = 1 .. nops(timePoints))]); if not 
SysInput = [] then inputTimes := SysInput:-getAllValuesOfVariable(time); 
inputValues := SysInput:-getAllValuesOfVariable(NameInputParameter); 
entriesSmallerLastTimePoint := select(`<`,convert(inputTimes,list),endTime); 
positionJustSmallerThanEnd := ListTools[BinarySearch](convert(inputTimes,list)
,entriesSmallerLastTimePoint[-1]); SetActInputPosition(
positionJustSmallerThanEnd+1); inputTable := table([seq(evalf(inputTimes[i]) =
inputValues[i],i = 1 .. positionJustSmallerThanEnd)]); integrationTime := [
indices(timeTable,nolist)]; for actInputIndex in indices(inputTable,nolist) 
while true do evalbList := map(evalb,[seq(actInputIndex = integrationTime[i],i
= 1 .. nops(integrationTime))]); if member(true,evalbList,'indexTrue') then 
timeTable[integrationTime[indexTrue]] := 3; actInputValue := inputTable[
actInputIndex]; inputTable[actInputIndex] := 'inputTable[%]'; inputTable[
integrationTime[indexTrue]] := actInputValue else timeTable[actInputIndex] :=
2 end if end do; allTimePoints := sort([indices(timeTable,nolist)]) else 
allTimePoints := timePoints end if; integrationCounter := 0; for k to nops(
allTimePoints) do actTOUT := allTimePoints[k]; setEndTime(actTOUT); actResult
:= integrate(printInfo); statusInt := actResult[2]; if not statusInt = 3 then
error "Integrating to TOUT=%1 was not succesfull. DDASPKSen terminated with ID\
ID = %2", actTOUT, statusInt end if; if timeTable[actTOUT] = 2 or timeTable[
actTOUT] = 3 then setSpecificParameters([NameInputParameter = inputTable[
actTOUT]]); calculateYPRIME(); numberOfVariables := nops([op(Sys[AlgVars]), op
(Sys[DynVars])]); for i from numberOfVariables+1 to NEQ do y_vec[i] := 0 end 
do; INFO_vec[1] := 0 end if; if not timeTable[actTOUT] = 2 then 
integrationCounter := integrationCounter+1; resultTable[time][
integrationCounter+1] := actTOUT; actData := actResult[1]; for item in 
varNames while true do resultTable[item][integrationCounter+1] := actData[item
] end do end if end do; otherData := [('NumPoints') = nrTimePoints+1, 
Parameters = Sys["Parameters"], DynVars = Sys["DynVars"], SenPars = Sys[
"SenPars"], AlgVars = Sys["AlgVars"], ExplicitAEs = Sys["ExplicitAEs"]]; curve
:= Aux:-CreateCurve(resultTable,otherData); return eval(curve) end proc; 
runInt := proc () local NumVarsInDAESys, Pars, ValuesExplicitVars, OtherData,
SenNames, Curve; SenNames := map(op,[seq([seq(convert(cat(d,convert(Sys[
"DynVars"][i2],string),d,convert(lhs(Sys["SenPars"][i1]),string)),name),i2 = 1
.. nops(Sys["DynVars"])), seq(convert(cat(d,convert(Sys["AlgVars"][i2],string)
,d,convert(lhs(Sys["SenPars"][i1]),string)),name),i2 = 1 .. nops(Sys["AlgVars"
]))],i1 = 1 .. nops(Sys["SenPars"]))]); VarNames := [op(Sys["DynVars"]), op(
Sys["AlgVars"]), op(SenNames), op(map(lhs,Sys["ExplicitAEs"]))]; NumVars := 
nops(VarNames); data := table(); for item in VarNames while true do data[item]
:= array(1 .. MaxNumSteps) end do; if member(t,VarNames) then error 
"name t must not be a variable names" end if; data[t] := array(1 .. 
MaxNumSteps); for i1 to MaxNumSteps do if `mod`(i1,100) = 0 then printf(
"done with %d out of %d steps\n",i1,MaxNumSteps) end if; DdaspkStatus := 
doIntStep(); ValuesExplicitVars := map(evalf,map(rhs,getExplicitAEs())); 
ValuesExplicitVars := subs(t = getT(),ValuesExplicitVars); CurrentPoint := [op
(getY()), op(ValuesExplicitVars)]; for i2 to nops(VarNames) do item := 
VarNames[i2]; data[item][i1] := CurrentPoint[i2] end do; data[t][i1] := getT()
; if not DdaspkStatus = 1 then break end if end do; OtherData := [('NumPoints'
) = i1-1]; OtherData := [op(OtherData), DDASPKstatus = DdaspkStatus]; 
OtherData := [op(OtherData), Parameters = [op(getParameters()), op(
getSenParameters())]]; Curve := Aux:-CreateCurve(data,OtherData); return Curve
end proc; runIntWithDynConstraint := proc (DynCon::term, ReqEps::EvalsToFloat,
ReqKeyword::name) local NumVarsInDAESys, Pars, ValuesExplicitVars, OtherData,
SenNames, Curve, CurrentSubsList, ValueOfDynCon, OldValueOfDynCon, Eps, 
NumOfCriticalPoint, MaxDeviation, TableCriticalPoints, CurrentMaxPoint, 
NumOfEntriesGreaterThanCurrentMin; OldValueOfDynCon := .1e21; 
NumOfCriticalPoint := 0; MaxDeviation := false; 
NumOfEntriesGreaterThanCurrentMin := 0; if 2 <= nargs then Eps := ReqEps else
Eps := .1e-2 end if; if nargs = 3 then if not ReqKeyword = ('maximum') then 
printf("WARNING: Keyword %1 must be 'maximum', will be ignored",ReqKeyword) 
else MaxDeviation := true; TableCriticalPoints := table() end if end if; 
SenNames := [seq(seq(cat(d,Sys["DynVars"][i1],d,lhs(Sys["SenPars"][i2])),i2 =
1 .. nops(Sys["SenPars"])),i1 = 1 .. nops(Sys["DynVars"]))]; VarNames := [op(
Sys["DynVars"]), op(Sys["AlgVars"]), op(SenNames), op(map(lhs,Sys[
"ExplicitAEs"]))]; NumVars := nops(VarNames); data := table(); for item in 
VarNames while true do data[item] := array(1 .. MaxNumSteps) end do; if member
(t,VarNames) then error "name t must not be a variable names" end if; data[t]
:= array(1 .. MaxNumSteps); for i1 to MaxNumSteps do if `mod`(i1,100) = 0 then
printf("done with %d out of %d steps\n",i1,MaxNumSteps) end if; DdaspkStatus 
:= doIntStep(); ValuesExplicitVars := map(evalf,map(rhs,getExplicitAEs())); 
ValuesExplicitVars := subs(t = getT(),ValuesExplicitVars); CurrentPoint := [op
(getY()), op(ValuesExplicitVars)]; for i2 to nops(VarNames) do item := 
VarNames[i2]; data[item][i1] := CurrentPoint[i2] end do; data[t][i1] := getT()
; if not MaxDeviation and NumOfCriticalPoint = 0 then CurrentSubsList := [t =
getT(), seq(VarNames[i2] = CurrentPoint[i2],i2 = 1 .. nops(VarNames))]; 
ValueOfDynCon := evalf(subs(CurrentSubsList,DynCon)); if ValueOfDynCon < -Eps
then printf("The dynamic constraint %f is smaller than the tolerance %f \n",
ValueOfDynCon,Eps); NumOfCriticalPoint := NumOfCriticalPoint+1 end if elif 
MaxDeviation then CurrentSubsList := [t = getT(), seq(VarNames[i2] = 
CurrentPoint[i2],i2 = 1 .. nops(VarNames))]; ValueOfDynCon := evalf(subs(
CurrentSubsList,DynCon)); if OldValueOfDynCon < ValueOfDynCon then if 
NumOfEntriesGreaterThanCurrentMin = 0 then printf(
"The dynamic constraint %f is smaller than the tolerance %f \n",
OldValueOfDynCon,Eps); NumOfCriticalPoint := NumOfCriticalPoint+1; 
TableCriticalPoints[NumOfCriticalPoint] := CurrentMaxPoint end if; 
NumOfEntriesGreaterThanCurrentMin := NumOfEntriesGreaterThanCurrentMin+1; 
OldValueOfDynCon := ValueOfDynCon end if; if ValueOfDynCon < -Eps and 
ValueOfDynCon < OldValueOfDynCon then OldValueOfDynCon := ValueOfDynCon; 
CurrentMaxPoint := CurrentSubsList; NumOfEntriesGreaterThanCurrentMin := 0 end
if end if; if not DdaspkStatus = 1 then if NumOfCriticalPoint = 0 then 
CurrentSubsList := []; TableCriticalPoints := table() end if; break end if end
do; OtherData := [('NumPoints') = i1-1]; OtherData := [op(OtherData), 
DDASPKstatus = DdaspkStatus]; print("RunInt finishes with DDASPK status ",
DdaspkStatus); OtherData := [op(OtherData), Parameters = GetParameters()]; 
Curve := Aux:-CreateCurve(data,OtherData); if not MaxDeviation then return 
Curve, CurrentSubsList else return Curve, TableCriticalPoints end if end proc;
setEndTime := proc (EndTime::EvalsToFloat) tout_vec := eval(EndTime); return 
end proc; setInitialTime := proc (InitialTime::EvalsToFloat) t_vec := eval(
InitialTime); return end proc; setMaxNumSteps := proc (RequestedMaxNumSteps::
posint) MaxNumSteps := RequestedMaxNumSteps; return end proc; setParameters :=
proc (NewPars::list(name = EvalsToFloat)) local AssignedPars, MissingAssigns,
ParNamesOfInstance, ListOfParValues; AssignedPars := map(lhs,NewPars); 
ParNamesOfInstance := map(lhs,Sys["Parameters"]); MissingAssigns := `minus`(
convert(ParNamesOfInstance,set),convert(AssignedPars,set)); if not 
MissingAssigns = {} then error cat(
"all parameters must be assigned by first argument, parameters %1 are missing"
,MissingAssigns) end if; ListOfParValues := subs(NewPars,ParNamesOfInstance);
setPars(ListOfParValues); return end proc; setPars := proc (NewPars::list(
EvalsToFloat)) local ListOfParamNames; if not nops(NewPars) = NumParsInDAEsys
then error "DAE system has %1 parameters", NumParsInDAEsys end if; rpar_vec :=
Vector(1 .. NumParsInDAEsys,map(evalf,NewPars),datatype = float[8]); 
ListOfParamNames := map(lhs,Sys["Parameters"]); Sys["Parameters"] := [seq(
ListOfParamNames[i1] = NewPars[i1],i1 = 1 .. nops(ListOfParamNames))]; return
end proc; setSenParameters := proc (NewPars::list(name = EvalsToFloat)) local
AssignedPars, MissingAssigns, SenParNamesOfInstance, ListOfParValues; 
AssignedPars := map(lhs,NewPars); SenParNamesOfInstance := map(lhs,Sys[
"SenPars"]); MissingAssigns := `minus`(convert(SenParNamesOfInstance,set),
convert(AssignedPars,set)); if not MissingAssigns = {} then error cat("all sen\
sitivity parameters must be assigned by first argument, sensitivity parameters\
 %1 are missing",MissingAssigns) end if; ListOfParValues := subs(NewPars,
SenParNamesOfInstance); setSenPars(ListOfParValues); return end proc; 
setSenPars := proc (NewPars::list(EvalsToFloat)) local ListOfParamNames; if 
not nops(NewPars) = NumSenParsInDAESys then error 
"DAE system has %1 sensitivity parameters", NumSenParsInDAESys end if; 
senpar_vec := Vector(1 .. NumSenParsInDAESys,map(evalf,NewPars),datatype = 
float[8]); ListOfParamNames := map(lhs,Sys["SenPars"]); Sys["SenPars"] := [seq
(ListOfParamNames[i1] = NewPars[i1],i1 = 1 .. nops(ListOfParamNames))]; return
end proc; setVariables := proc (NewVars::list(name = EvalsToFloat)) local 
AssignedVars, MissingAssigns, VarNamesOfInstance, ListOfVarValues; 
AssignedVars := map(lhs,NewVars); VarNamesOfInstance := [op(Sys["DynVars"]), 
op(Sys["AlgVars"])]; MissingAssigns := `minus`(convert(VarNamesOfInstance,set)
,convert(AssignedVars,set)); if not MissingAssigns = {} then error cat(
"all variables must be assigned by first argument, variables %1 are missing",
MissingAssigns) end if; ListOfVarValues := subs(NewVars,VarNamesOfInstance); 
SetVars(ListOfVarValues); return end proc; setVars := proc (NewVars::list(
EvalsToFloat)) if not nops(NewVars) = nops(Sys["DynVars"])+nops(Sys["AlgVars"]
) then error "DAE system has %1 Variables", nops(Sys[DynVars])+nops(Sys[
AlgVars]) end if; y_vec := Vector(1 .. NEQ,NewVars,datatype = float[8]); 
return end proc; createCFunction := proc () local i, fh, listT_vec, listY_vec,
listYprime_vec, listTout_vec, listINFO_vec, listRtol_vec, listAtol_vec, 
listRWORK_vec, listIWORK_vec, listRpar_vec, listIpar_vec, listSenpar_vec, 
compileCommand, mpiCompileCommand, PATH, status, NumVarsInDAESys, Pars, 
ValuesExplicitVars, OtherData, SenNames, Curve, z; currentdir("simModel"); 
listT_vec := convert(eval(t_vec),list); listY_vec := convert(eval(y_vec),list)
; listYprime_vec := convert(eval(yprime_vec),list); listTout_vec := convert(
eval(tout_vec),list); listINFO_vec := convert(eval(INFO_vec),list); 
listRtol_vec := convert(eval(rtol_vec),list); listAtol_vec := convert(eval(
atol_vec),list); listRWORK_vec := convert(eval(RWORK_vec),list); listIWORK_vec
:= convert(eval(IWORK_vec),list); listRpar_vec := convert(eval(rpar_vec),list)
; listIpar_vec := convert(eval(ipar_vec),list); listSenpar_vec := convert(eval
(senpar_vec),list); status := ssystem("rm -fr ./cwrapperDDASPK.c"); if status[
1] <> 0 then error " in \".../DDASPKSen/source/SubsIntoSimpleWrapper/SimpleWra\
pper/createCFunction.txt\" command:\"rm -fr ./ParamEstWrapper.c\"" end if; fh
:= fopen("cwrapperDDASPK.c",WRITE); fprintf(fh,"#include <stdio.h>\n\n"); 
fprintf(fh,"extern \"C\" void cwrapper_ddaspk(\n"); fprintf(fh,"  int, \n"); 
fprintf(fh,"  double *, \n"); fprintf(fh,"  double *, \n"); fprintf(fh,
"  double *, \n"); fprintf(fh,"  double , \n"); fprintf(fh,"  int *, \n"); 
fprintf(fh,"  double *, \n"); fprintf(fh,"  double *, \n"); fprintf(fh,
"  int *, \n"); fprintf(fh,"  double *, \n"); fprintf(fh,"  int, \n"); fprintf
(fh,"  int *, \n"); fprintf(fh,"  int, \n"); fprintf(fh,"  double *, \n"); 
fprintf(fh,"  int, \n"); fprintf(fh,"  double *\n"); fprintf(fh,");\n\n"); 
fprintf(fh,
"void cwrapperDDASPK(double* rpar_vec,double* data_vec, double* y_vec){\n\n");
fprintf(fh,"int    NEQ = %d;\n",NEQ); fprintf(fh,"double t_vec[1]={0};\n"); 
fprintf(fh,cat("double tout_vec[%d]={",seq("%1.20e, ",i = 1 .. nops(
listTout_vec)-1),"%1.20e};\n"),nops(listTout_vec),op(listTout_vec)); fprintf(
fh,cat("int INFO_vec[%d]={",seq("%1.20e, ",i = 1 .. nops(listINFO_vec)-1),
"%1.20e};\n"),nops(listINFO_vec),op(listINFO_vec)); fprintf(fh,
"double rtol_vec[%d]={%de%d};\n",nops(listRtol_vec)-1,op(listRtol_vec)); 
fprintf(fh,"double atol_vec[%d]={%de%d};\n",nops(listAtol_vec)-1,op(
listAtol_vec)); fprintf(fh,"int    idid=0;\n"); fprintf(fh,cat(
"double RWORK_vec[%d]={",seq("%1.20e, ",i = 1 .. nops(listRWORK_vec)-1),
"%1.20e};\n"),nops(listRWORK_vec),op(listRWORK_vec)); fprintf(fh,
"int    LRW = %d;\n",LRW); fprintf(fh,cat("int IWORK_vec[%d]={",seq("%1.20e, "
,i = 1 .. nops(listIWORK_vec)-1),"%1.20e};\n"),nops(listIWORK_vec),op(
listIWORK_vec)); fprintf(fh,"int    LIW = %d;\n",LIW); fprintf(fh,cat(
"int ipar_vec[%d]={",seq("%1.20e, ",i = 1 .. nops(listIpar_vec)-1),"%1.20e};\n\
"),nops(listIpar_vec),op(listIpar_vec)); fprintf(fh,cat("double yprime_vec[",
NEQ,"];")); if 0 < nops(listSenpar_vec) then fprintf(fh,cat(
"\ndouble senpar_vec[%d]={",seq("%1.20e, ",i = 1 .. nops(listSenpar_vec)-1),
"%1.20e};\n"),nops(listSenpar_vec),op(listSenpar_vec)) else fprintf(fh,
"\n\ndouble *senpar_vec;") end if; fprintf(fh,"\ndouble tmpData_vec[%d*%d];",
MaxNumSteps,nops(listY_vec)+1); fprintf(fh,"\n\n tmpData_vec[0]=t_vec[0];"); 
fprintf(fh,"\n\nfor(int i=0; i<%d; i++){\n",nops(listY_vec)); fprintf(fh,
"tmpData_vec[1+i]=y_vec[i];\n"); fprintf(fh,"}\n"); fprintf(fh,"\nint steps;\n\
"); fprintf(fh,cat("\nfor(steps=1; steps<=",MaxNumSteps,";steps++){\n")); 
fprintf(fh,"\n\n \n\ncwrapper_ddaspk(NEQ,t_vec,y_vec,yprime_vec,tout_vec[0],IN\
FO_vec,rtol_vec,atol_vec,&idid,RWORK_vec,LRW,IWORK_vec,LIW,rpar_vec,ipar_vec[0\
],senpar_vec);"); fprintf(fh,"\n\ntmpData_vec[steps*%d]=t_vec[0];",nops(
listY_vec)+1); fprintf(fh,"\n\n for(int j=0; j<%d; j++){",nops(listY_vec)); 
fprintf(fh,"\n  tmpData_vec[steps*%d+j+1]=y_vec[j];",nops(listY_vec)+1); 
fprintf(fh,"\n }\n"); z := 0; for z to nops(VarNames)-nops(listY_vec) do 
fprintf(fh,StringTools[Remove](isApostrophe(z),convert(cat("\ndata[j][i]=",
convert(rhs(subs(VarSeq,ParSeq,Sys["ExplicitAEs"][z])),string),";"),string)));
fprintf(fh,"\nprintf(\"data[%%i][%%i]=%%f   \",j,i,data[j][i]); "); fprintf(fh
,"\nj++;") end do; fprintf(fh,"\n\n"); fprintf(fh,
"\n if(idid==-33){\nINFO_vec[0]=0;\nINFO_vec[2]=1;\n"); fprintf(fh,
"\n  printf(\"\\n\\n ERROR idid(-33) ERROR\");"); fprintf(fh,"\n \n\n \n\ncwra\
pper_ddaspk(NEQ,t_vec,y_vec,yprime_vec,tout_vec[0],INFO_vec,rtol_vec,atol_vec,\
&idid,RWORK_vec,LRW,IWORK_vec,LIW,rpar_vec,ipar_vec[0],senpar_vec);"); fprintf
(fh,"\n  printf(\"\\n\\nnew idid= %%d\",idid);\n }"); fprintf(fh,
"\n\n if(idid!=1) break;"); fprintf(fh,"\n}"); fprintf(fh,"\nsteps++;"); 
fprintf(fh,"\nfor(int i=0;i<%d;i++){",nops(listY_vec)+1); fprintf(fh,
"\n\n for(int j=0;j<steps;j++){"); fprintf(fh,
"\n  data_vec[i*steps+j+1] = tmpData_vec[j*%d+i];",nops(listY_vec)+1); fprintf
(fh,"\n }"); fprintf(fh,"\n}\n"); fprintf(fh,"data_vec[0]=steps;\n"); fprintf(
fh,"}\n"); fclose(fh); fh := fopen("settings.txt",WRITE); fprintf(fh,"%f\n",
nops(listY_vec)+nops(listRpar_vec)); fprintf(fh,"%f\n",nops(listY_vec)); 
fprintf(fh,"%f\n",MaxNumSteps); fclose(fh); status := ssystem(cat("rm -fr ",
_EnvLibDir,"cwrap_ddaspk?.so")); printf("%a",currentdir()); printf("%a",cat(
"cp *.so ",_EnvLibDir,"/")); status := ssystem(cat("cp cwrap_ddaspk",convert(
InstanceCounter,string),".so ",_EnvLibDir,"/")); if status[1] <> 0 then error
" in \".../DDASPKSen/source/SubsIntoSimpleWrapper/SimpleWrapper/createCFunctio\
n.txt\"command:\"cp cwrap_ddaspk\",convert(InstanceCounter, string),\".so \""
end if; currentdir(".."); RETURN() end proc; isApostrophe := proc (s) evalb(s
= "`") end proc; calculateIntResults := proc (timePoints::list(numeric)) local
intiY, initTime, initPoint, timeIsOrdered, senNames, varNames, resultTable, 
results, nrTimePoints, actData, actName, actTOUT, actResult, statusInt, 
lastPosition, otherData, curve, valuesExplicitVars, i, inputTimes, inputValues
, intAndInputTimes, k, initY, positionJustSmallerThanEnd, 
entriesSmallerLastTimePoint, endTime, sysVarNames, timeTable, inputTable, 
integrationTime, actInputIndex, evalbList, actInputValue, allTimePoints, 
integrationCounter; nrTimePoints := nops(timePoints); endTime := timePoints[-1
]; initY := getY(); initTime := getT(); valuesExplicitVars := map(evalf,map(
rhs,getExplicitAEs())); CheckInputForCalculateIntResults(timePoints,initTime);
if not type(initY,list(numeric)) then error 
"Initial values for variables have not been set. Call SetVariables()." end if;
calculateYPRIME(); INFO_vec[1] := 0; senNames := map(op,[seq([seq(convert(cat(
d,convert(Sys["DynVars"][i2],string),d,convert(lhs(Sys["SenPars"][i1]),string)
),name),i2 = 1 .. nops(Sys["DynVars"])), seq(convert(cat(d,convert(Sys[
"AlgVars"][i2],string),d,convert(lhs(Sys["SenPars"][i1]),string)),name),i2 = 1
.. nops(Sys["AlgVars"]))],i1 = 1 .. nops(Sys["SenPars"]))]); sysVarNames := [
op(Sys["DynVars"]), op(Sys["AlgVars"])]; varNames := [op(sysVarNames), op(
senNames), op(map(lhs,Sys["ExplicitAEs"]))]; resultTable := table(); 
resultTable[time] := array(1 .. nrTimePoints+1); resultTable[time][1] := 
initTime; initPoint := [op(initY), op(valuesExplicitVars)]; for i to nops(
varNames) do actName := varNames[i]; resultTable[actName] := array(1 .. 
nrTimePoints+1); resultTable[actName][1] := initPoint[i] end do; timeTable :=
table([seq(evalf(timePoints[i]) = 1,i = 1 .. nops(timePoints))]); if not 
SysInput = [] then inputTimes := SysInput:-getAllValuesOfVariable(time); 
inputValues := SysInput:-getAllValuesOfVariable(NameInputParameter); 
entriesSmallerLastTimePoint := select(`<`,convert(inputTimes,list),endTime); 
positionJustSmallerThanEnd := ListTools[BinarySearch](convert(inputTimes,list)
,entriesSmallerLastTimePoint[-1]); SetActInputPosition(
positionJustSmallerThanEnd+1); inputTable := table([seq(evalf(inputTimes[i]) =
inputValues[i],i = 1 .. positionJustSmallerThanEnd)]); integrationTime := [
indices(timeTable,nolist)]; for actInputIndex in indices(inputTable,nolist) 
while true do evalbList := map(evalb,[seq(actInputIndex = integrationTime[i],i
= 1 .. nops(integrationTime))]); if member(true,evalbList,'indexTrue') then 
timeTable[integrationTime[indexTrue]] := 3; actInputValue := inputTable[
actInputIndex]; inputTable[actInputIndex] := 'inputTable[%]'; inputTable[
integrationTime[indexTrue]] := actInputValue else timeTable[actInputIndex] :=
2 end if end do; allTimePoints := sort([indices(timeTable,nolist)]) else 
allTimePoints := timePoints end if; integrationCounter := 0; for k to nops(
allTimePoints) do actTOUT := allTimePoints[k]; setEndTime(actTOUT); actResult
:= runInt(); statusInt := rhs(actResult:-getNonTableData()[2]); if not 
statusInt = 3 then error "Integrating to TOUT=%1 was not succesfull. DDASPKSen\
 terminated with IDID = %2", actTOUT, statusInt end if; if timeTable[actTOUT]
= 2 or timeTable[actTOUT] = 3 then setSpecificParameters([NameInputParameter =
inputTable[actTOUT]]); calculateYPRIME(); INFO_vec[1] := 0 end if; if not 
timeTable[actTOUT] = 2 then integrationCounter := integrationCounter+1; 
lastPosition := nops(op(3,eval(actResult:-getAllValuesOfVariable(t)))); 
resultTable[time][integrationCounter+1] := actTOUT; actData := actResult:-
getData(); for item in varNames while true do resultTable[item][
integrationCounter+1] := actData[item][lastPosition] end do end if end do; 
otherData := [('NumPoints') = nrTimePoints+1, Parameters = Sys["Parameters"],
DynVars = Sys["DynVars"], SenPars = Sys["SenPars"], AlgVars = Sys["AlgVars"],
ExplicitAEs = Sys["ExplicitAEs"]]; curve := Aux:-CreateCurve(resultTable,
otherData); return eval(curve) end proc; setInput := proc (inputCurve::Curve)
local inputVarNames, actParName, timePoints, timeIsOrdered, actTime, allPars;
inputVarNames := inputCurve:-getVarNames(); if not (ActInputPosition = 0 and 
SysInput = [] and NameInputParameter = 0) then warning(
"Input has alread been set. Old settings are overwritten by new ones.") end if
; if not `in`(time,inputVarNames) then error 
"The inputCurve does not contain time as a variable" elif not nops(
inputVarNames) = 2 then error "The inputCurve has to contain two variables: ti\
me, and the parameter describing the input" end if; actParName := 
inputVarNames[-ListTools[Search](time,inputVarNames)]; allPars := map(lhs,
aDAESys["Parameters"]); if not `in`(actParName,allPars) then error "The input \
parameter \"%1\" has to be a parameter from the list of parameters specified i\
n this DAE System \"%2\"", actParName, allPars end if; timePoints := convert(
inputCurve:-getAllValuesOfVariable(time),list); timeIsOrdered := ListTools[
Sorted](timePoints,verify,less_than); actTime := getT(); if not timeIsOrdered
then error 
"Input time points have to be strictly ascending. You have inserted: %1", 
timePoints elif not actTime <= timePoints[1] then error "The first time point \
of the time series has to be >= the integration timepoint while SetInput() is \
envoked" end if; SysInput := inputCurve; NameInputParameter := actParName; 
ActInputPosition := 1 end proc; getInputCurve := proc () if SysInput = [] then
error "Input is not set. This can be done with SetInput()" else return eval(
SysInput) end if end proc; getNameInputParameter := proc () if 
NameInputParameter = 0 then error 
"Input parameter is not set. This can be done with SetInput()" else return 
eval(NameInputParameter) end if end proc; getActInputPosition := proc () if 
ActInputPosition = 0 then Warning(
"Input has not yet been set. This can be done with setInput().") else return 
ActInputPosition end if end proc; setSpecificParameters := proc (NewPars::list
(name = EvalsToFloat)) local AssignedPars, ParNamesOfInstance, 
ParAssignmentsInstance, ListOfParValues, actAsssignment, actIndex, actParName;
AssignedPars := map(lhs,NewPars); ParAssignmentsInstance := Sys["Parameters"];
ParNamesOfInstance := map(lhs,ParAssignmentsInstance); if not `subset`(convert
(AssignedPars,set),convert(ParNamesOfInstance,set)) then error cat("parameter \
names in first argument must be parameter names from the DAESys: %1.",
ParNamesOfInstance) end if; for actParName in AssignedPars while true do 
actIndex := ListTools[Search](actParName,ParNamesOfInstance); 
ParAssignmentsInstance[actIndex] := NewPars[actIndex] end do; ListOfParValues
:= map(rhs,ParAssignmentsInstance); SetPars(ListOfParValues); return end proc;
setActInputPosition := proc (aPosition::integer) if ActInputPosition = 0 then
Warning("Input has not yet been set.Therefore the ActInputPosition cannot be r\
eset. Set the input with SetInput().") else ActInputPosition := aPosition end
if end proc end module; RETURN(eval(SimpleWrapper)) end proc; Templates := 
module () export jacAdiforTemplate, jacTemplate, resTemplate; 
jacAdiforTemplate := proc (t::numeric, y, yprime, pd, cj::numeric, par, ipar,
senpar, ijac::integer) local UnitNumber, NumVarsTotal, NumEqnsTotal, ires, i1,
i2, ySeedMatrix, yprimeSeedMatrix, gdelta, delta, ldg_y, ldg_yprime, ldg_delta
, g_p, ehsupPar; declare(y = array(1 .. NumODEsToBeSubs+NumAEsToBeSubs,numeric
),yprime = array(1 .. NumODEsToBeSubs+NumAEsToBeSubs,numeric),pd = array(1 ..
NumODEsToBeSubs+NumAEsToBeSubs,1 .. NumDynVarsToBeSubs+NumAlgVarsToBeSubs,
numeric),par = array(1 .. NumParsToBeSubs,numeric),ipar = array(1 .. 1,numeric
),senpar = array(1 .. NumSenParsToBeSubs,numeric),ires = integer,ySeedMatrix =
array(1 .. 2*NumDynVarsToBeSubs+2*NumAlgVarsToBeSubs,1 .. NumDynVarsToBeSubs+
NumAlgVarsToBeSubs,numeric),yprimeSeedMatrix = array(1 .. 2*NumDynVarsToBeSubs
+2*NumAlgVarsToBeSubs,1 .. NumDynVarsToBeSubs+NumAlgVarsToBeSubs,numeric),
delta = array(1 .. NumODEsToBeSubs+NumAEsToBeSubs,numeric),gdelta = array(1 ..
2*NumDynVarsToBeSubs+2*NumAlgVarsToBeSubs,1 .. NumODEsToBeSubs+NumAEsToBeSubs,
numeric)); NumVarsTotal := NumDynVarsToBeSubs+NumAlgVarsToBeSubs; NumEqnsTotal
:= NumODEsToBeSubs+NumAEsToBeSubs; UnitNumber := UnitNumberToBeSubs; ehsupPar
:= -1; for i1 to NumVarsTotal do for i2 to NumVarsTotal do ySeedMatrix[i1,i2]
:= 0.; ySeedMatrix[NumVarsTotal+i1,i2] := 0.; yprimeSeedMatrix[i1,i2] := 0.; 
yprimeSeedMatrix[NumVarsTotal+i1,i2] := 0. end do; ySeedMatrix[i1,i1] := 1.0;
yprimeSeedMatrix[NumVarsTotal+i1,i1] := 1.0 end do; ldg_y := 2*NumVarsTotal; 
ldg_yprime := 2*NumVarsTotal; ldg_delta := 2*NumVarsTotal; g_p := 2*
NumVarsTotal; NameOfDerivateRoutineToBeSubs(g_p,t,y,ySeedMatrix,ldg_y,yprime,
yprimeSeedMatrix,ldg_yprime,cj,delta,gdelta,ldg_delta,ires,par,ipar,senpar); 
for i1 to NumEqnsTotal do for i2 to NumVarsTotal do pd[i1,i2] := gdelta[i2,i1]
+cj*gdelta[i2+NumVarsTotal,i1] end do end do; RETURN() end proc; jacTemplate 
:= proc (t::numeric, y, yprime, pd, cj::numeric, par, ipar, senpar, ijac::
integer) declare(y = array(1 .. NumODEsToBeSubs+NumAEsToBeSubs,numeric),yprime
= array(1 .. NumODEsToBeSubs+NumAEsToBeSubs,numeric),pd = array(1 .. 
NumODEsToBeSubs+NumAEsToBeSubs,1 .. NumDynVarsToBeSubs+NumAlgVarsToBeSubs,
numeric),par = array(1 .. NumParsToBeSubs,numeric),ipar = array(1 .. 1,numeric
),senpar = array(1 .. NumSenParsToBeSubs,numeric)); pd := array(1 .. 
NumODEsToBeSubs+NumAEsToBeSubs,1 .. NumDynVarsToBeSubs+NumAlgVarsToBeSubs,
JacToBeSubs); RETURN() end proc; resTemplate := proc (t::numeric, y, yprime, 
cj::numeric, delta, ires::integer, par, ipar, senpar) local z; declare(y = 
array(1 .. NumODEsToBeSubs+NumAEsToBeSubs,numeric),yprime = array(1 .. 
NumODEsToBeSubs+NumAEsToBeSubs,numeric),delta = array(1 .. NumODEsToBeSubs+
NumAEsToBeSubs,numeric),par = array(1 .. NumParsToBeSubs,numeric),ipar = array
(1 .. 1,integer),senpar = array(1 .. NumSenParsToBeSubs,numeric),z = array(1 
.. NumExplicitAEsToBeSubs+1,numeric)); z := ExplicitAEsToBeSubs; delta := 
ResToBeSubs; RETURN() end proc; end module; InstanceCounter := 0; 
appendSenParsToDAESys := proc (aDAESys::DAESys, SensitivityPars::{name, list(
name), list(name = EvalsToFloat), name = EvalsToFloat}) local SenParsList, 
SetOfIndices, NewDAESys, SetOfCommonNames, i1; NewDAESys := copy(aDAESys); if
type(SensitivityPars,name = EvalsToFloat) or type(SensitivityPars,name) then 
SenParsList := convert(SensitivityPars,list) else SenParsList := 
SensitivityPars end if; if type(SenParsList,list(name)) then SenParsList := 
Aux:-ListOperations:-subsToCreateSubsList(NewDAESys["Parameters"],SenParsList)
end if; SetOfIndices := {indices(NewDAESys)}; if member(["SenPars"],
SetOfIndices) then SetOfCommonNames := `intersect`(map(lhs,convert(SenParsList
,set)),map(lhs,convert(NewDAESys["SenPars"],set))); if not SetOfCommonNames =
{} then printf("Warning: %q already in [SenPars]\n",op(SetOfCommonNames)); 
SenParsList := Aux:-ListOperations:-removeItemFromList([op(SetOfCommonNames)],
SenParsList) end if else NewDAESys["SenPars"] := [] end if; for i1 to nops(
SenParsList) do if Aux:-ListOperations:-isLHSin(lhs(SenParsList[i1]),NewDAESys
["Parameters"]) then NewDAESys["SenPars"] := [op(NewDAESys["SenPars"]), 
SenParsList[i1]]; NewDAESys["Parameters"] := Aux:-ListOperations:-
removeItemFromList(lhs(SenParsList[i1]),NewDAESys["Parameters"]) else printf(
"WARNING: %s is not a parameter of the DAE system, \n",lhs(SenParsList[i1]));
printf("it will be ignored.") end if end do; return eval(NewDAESys) end proc;
init := proc () InstanceCounter := 0 end proc; CreateInstance := proc (
ReqDAESys::DAESys, WorkingDir::string, SenParameters::{name, list(name), list(
name = EvalsToFloat), name = EvalsToFloat}) local ddaspk, NEQtoBeSubs, 
LIWtoBeSubs, LRWtoBeSubs, NY, base, SimpleWrapper, INFO5toBeSubs, aDAESys, 
DAEsystem, NumParsInDAEsys, NumSenParsInDAESys, NewDAESys, NamesDynVars, 
NamesAlgVars, NamesExplicitAEs, InterfaceKeyword, NamesPars, NamesSenPars, 
KnownOpts, UnknownOpts, ReqOpts, NameOfJac, NameOfRes, SetEndOfFileName, Jac,
Res; InstanceCounter := InstanceCounter+1; aDAESys := appendSenParsToDAESys(
ReqDAESys,SenParameters); KnownOpts := ['adifor', 'NoInterface', 'SetFileName'
]; if 3 < nargs then ReqOpts := [args[4 .. -1]] else ReqOpts := [] end if; 
UnknownOpts := `minus`(convert(ReqOpts,set),convert(KnownOpts,set)); if not 
UnknownOpts = {} then WARNING(cat(
"in DDASPK:-CreateInstance, requested options %1 are ",
"not known. Known options are %2"),UnknownOpts,KnownOpts) end if; if member('
NoInterface',ReqOpts) then InterfaceKeyword := 'NoInterface' else 
InterfaceKeyword := 'BuildInterface' end if; if member('SetFileName',ReqOpts)
and InterfaceKeyword = ('NoInterface') then SetEndOfFileName := true else 
SetEndOfFileName := false end if; if 3 < nargs then if member('adifor',ReqOpts
) then UseAdifor := true end if else UseAdifor := false end if; if UseAdifor =
true then if not type(ADIFOR,`module`) then error 
"option adifor not available as module ADIFOR is not defined" end if end if; 
if member(["SenPars"],{indices(aDAESys)}) then EvalSen := true else EvalSen :=
false; aDAESys["SenPars"] := [] end if; if prepareDirForDDASPK(WorkingDir) = 
false then error "requested working dir %1 is already in use", WorkingDir end
if; NamesDynVars := {seq(op(1,aDAESys["DynVars"][i1]),i1 = 1 .. nops(aDAESys[
"DynVars"]))}; NamesAlgVars := {seq(op(1,aDAESys["AlgVars"][i1]),i1 = 1 .. 
nops(aDAESys["AlgVars"]))}; NamesExplicitAEs := {seq(op(1,map(rhs,aDAESys[
"ExplicitAEs"])[i1]),i1 = 1 .. nops(aDAESys["ExplicitAEs"]))}; NamesPars := {
seq(op(1,aDAESys["Parameters"][i1]),i1 = 1 .. nops(aDAESys["Parameters"]))}; 
NamesSenPars := {seq(op(1,aDAESys["SenPars"][i1]),i1 = 1 .. nops(aDAESys[
"SenPars"]))}; if NamesDynVars = {y} and (NamesAlgVars = {y} or NamesAlgVars =
{}) and (NamesExplicitAEs = {z} or NamesExplicitAEs = {}) and (NamesPars = {
par} or NamesPars = {} and NamesSenPars = {senpar} or NamesPars = {}) then 
WARNING("note that DDASPK:-CreateInstance() does not require system to be in s\
tandard notation"); NewDAESys := aDAESys else NewDAESys := Aux:-SystemClasses
:-subsStandardNotationIntoDAESysSen(aDAESys,y) end if; if SetEndOfFileName 
then Jac := convert(cat('jac',InstanceCounter),name); Res := convert(cat('res'
,InstanceCounter),name) else Jac := convert('jac',name); Res := convert('res',
name) end if; NameOfJac := convert(Jac,string); NameOfRes := convert(Res,
string); createExternalFunctions(eval(NewDAESys),Res,Jac); NumParsInDAEsys :=
nops(NewDAESys["Parameters"]); if EvalSen = true then NumSenParsInDAESys := 
nops(NewDAESys["SenPars"]) else NumSenParsInDAESys := 0 end if; ddaspk := 
createSharedObject(NameOfRes,Res,NameOfJac,Jac,WorkingDir,NewDAESys,
InterfaceKeyword); if InterfaceKeyword = ('BuildInterface') then DAEsystem :=
subs(map(rhs = lhs,NewDAESys["Substitutions"]),eval(NewDAESys)); DAEsystem[
"Substitutions"] := map(rhs = lhs,NewDAESys["Substitutions"]); SimpleWrapper 
:= SubsIntoSimpleWrapper(DAEsystem,ddaspk); save DAEsystem, cat(WorkingDir,
"/DAEsystem.m"); RETURN(eval(SimpleWrapper)) else RETURN(eval(aDAESys)) end if
end proc; RestoreInstance := proc (WorkingDir::string) local ddaspk, 
SharedObjectName, SharedObjectEnd, FileName, SimpleWrapper, OldDir, 
NumParsInDAEsys, NumSenParsInDAESys; SharedObjectEnd := ".so"; OldDir := 
currentdir(); currentdir(WorkingDir); FileName := "DAEsystem.m"; if not Aux:-
FileOperations:-fileExists(FileName) then currentdir(OldDir); error 
"in the requested dir %1, file named DAEsystem.m does not exist", WorkingDir 
end if; read "DAEsystem.m"; NumParsInDAEsys := nops(DAEsystem["Parameters"]);
NumSenParsInDAESys := nops(DAEsystem["Parameters"]); if 1 < nargs then if not
type(args[2],string) then error 
"optional third argument, SharedObjectName, must be a string" else 
SharedObjectName := args[2] end if else SharedObjectName := cat(
"cwrap_ddaspk1",SharedObjectEnd) end if; if not Aux:-FileOperations:-
fileExists(SharedObjectName) then currentdir(OldDir); error 
"in the requested dir %1, file named %2 does not exist", WorkingDir, 
SharedObjectName end if; if not 1 < nargs then warning(
"using default shared object cwrap_ddaspk1.so in %1",WorkingDir) end if; 
ddaspk := define_external('cwrapper_ddaspk',NEQ::integer[4],T::REF(float[8]),Y
::ARRAY(1 .. NEQ,float[8]),YPRIME::ARRAY(1 .. NEQ,float[8]),tout::float[8],
INFO::ARRAY(1 .. 30,integer[4]),RTOL::REF(float[8]),ATOL::REF(float[8]),IDID::
REF(integer[4]),RWORK::ARRAY(1 .. LENW,float[8]),LENW::integer[4],IWORK::ARRAY
(1 .. LENIW,integer[4]),LENIW::integer[4],RPAR::ARRAY(1 .. NumParsInDAEsys,
float[8]),IPAR::integer[4],SENPAR::ARRAY(1 .. NumSenParsInDAESys,float[8]),LIB
= SharedObjectName,WRAPPER); SimpleWrapper := SubsIntoSimpleWrapper(eval(
DAEsystem),ddaspk); currentdir(OldDir); unassign('DAEsystem'); RETURN(eval(
SimpleWrapper)) end proc; Tests := module () export all, simpleProblem1, 
simpleProblem1WithAdifor; all := proc () if simpleProblem1WithAdifor() = true
and simpleProblem1() = true then return true else return false end if end proc
; simpleProblem1 := proc () local g, xRange, yRange, anNLP, StandardNLP, 
NPSOLproc, Objf, OptPoint, ExpectedObjf, ExpectedOptPoint, Inform, aDAESys, 
StandardModel, DDASPKproc, OptValue, T, Y, TOUT, data, i, retVal, imax, 
result1, result2; ExpectedObjf := 1.00000; ExpectedOptPoint := [1.0, 1.0]; if
Aux:-FileOperations:-dirExists("tmp") then error 
"temporary directory ./tmp needed to run test already exists" else mkdir("tmp"
) end if; aDAESys := Aux:-SystemClasses:-newDAESys(); aDAESys["ODEs"] := [
`vm'` = (1-phis)*rhomf/rhom-vm*qoutoverqin-theta*vm*(lp+ltrm)*Q-2*theta*MWm/
rhom*f*kd*zi*If, `vs'` = phis*rhosf/rhos-vs*qoutoverqin, `zi'` = 1-theta*kd*zi
-zi*qoutoverqin, `y'` = rhof*Cpf/rho/Cp*(yf-yref)-rhof/rho*(y-yref)+negdeltaH/
rho/Cp/Tc*theta*vm*lp*rhom/MWm*Q-theta*alpha*B/rho/Cp*(y-1)]; aDAESys[
"Parameters"] := [theta = 200.4617, If = .316531e-1, phis = .6]; aDAESys[
"DynVars"] := [vm, vs, zi, y]; aDAESys["ExplicitAEs"] := [yf = 2/3, yref = 5/9
, Tc = 45, T = 45*y, negdeltaH = 21000, B = 1.4, alpha = 50/3, a3 = -.3495, a2
= -6.7530, a1 = -.4407, f = .8, R = 1.987, MWs = 74.10, MWm = 86.05, phim = 1-
phis, rhosf = 777.2814025, rhomf = 918.5720, rhomg = 892.0200, rhof = rhomf*
phim+rhosf*phis, rhos = 74120/(91.878+.116*T), rhop = 1211-.8496*T, rhom = 958\
.4-1.3276*T, rhofT = rhom*phim+rhos*phis, vp = 1-vm-vs, rho = rhom*vm+rhos*vs+
rhop*vp, kd = exp(ln(60)+34.99620124-30800/R/(T+273)), xt = (rhomf*phim-vm*
rhomg)/(rhomf*phim+MWm/MWs*phis*rhosf), lt = exp(a1*xt+a2*xt^2+a3*xt^3), ltrm
= 417.3948000*exp(-7569./R/(T+273.)), lp = 979.8*exp(-4869/R/(T+273))*15^(1/2)
, Cpf = (.470*rhomf*phim+.716*rhosf*phis)/(rhomf*phim+rhosf*phis), Cp = (.470*
rhom*vm+.716*rhos*vs+rhop*vp*(.321425+.955e-3*T))/(rhom*vm+rhos*vs+rhop*vp), Q
= 2^(1/2)*(f*kd*zi*If/lt)^(1/2), Rm = (lp+ltrm)*Q*rhom/MWm*vm+2*f*kd*zi*If, 
qoutoverqin = rhof/rhofT+theta*MWm*Rm*(1/rhop-1/rhom)]; aDAESys["AEs"] := [];
aDAESys["AlgVars"] := []; StandardModel := Aux:-SystemClasses:-
subsExplicitAEsIntoDAESys(aDAESys); DDASPKproc := DDASPKSen:-CreateInstance(
StandardModel,"./tmp",theta,'NoInterface'); OptValue := 1.25058770645712380; T
:= 0.; Y := [.1167175, .6434086, .9702816, 1.250588]; TOUT := 10.0; data := 
array(1 .. 1000,1 .. 5); DDASPKproc:-setVars(Y); DDASPKproc:-setInitialTime(T)
; DDASPKproc:-setEndTime(TOUT); DDASPKproc:-init(); for i to 999 do retVal :=
DDASPKproc:-doIntStep(); Y := DDASPKproc:-getY(); T := DDASPKproc:-getT(); 
data[i,1] := T; data[i,2] := Y[1]; data[i,3] := Y[2]; data[i,4] := Y[3]; data[
i,5] := Y[4]; if retVal = 1 then next else break end if end do; imax := i-1; 
for i to imax do if .1e-2 < abs(evalf(data[i,5])-OptValue) then error 
"test 1 using interface created by CreateInstance failed"; return false end if
end do; T := 0.; Y := [.5e-1, .6434086, .9702816, 1.250588]; TOUT := 10.0; 
data := array(1 .. 1000,1 .. 5); DDASPKproc:-setVars(Y); DDASPKproc:-
setInitialTime(T); DDASPKproc:-setEndTime(TOUT); DDASPKproc:-init(); for i to
999 do retVal := DDASPKproc:-doIntStep(); Y := DDASPKproc:-getY(); T := 
DDASPKproc:-getT(); data[i,1] := T; data[i,2] := Y[1]; data[i,3] := Y[2]; data
[i,4] := Y[3]; data[i,5] := Y[4]; if retVal = 1 then next else break end if 
end do; imax := i-1; if .1e-2 < abs(evalf(data[imax,5])-OptValue) then error 
"test 2 using interface created by CreateInstance failed"; return false end if
; try dlclose("cwrap_ddaspk1.so"); dlclose("mwrap_cwrapper_ddaspk.dll") catch:
NULL end try; DDASPKproc := 'DDASPKproc'; DDASPKproc := DDASPKSen:-
RestoreInstance("./tmp"); OptValue := 1.25058770645712380; T := 0.; Y := [.116\
7175, .6434086, .9702816, 1.250588]; TOUT := 10.0; data := array(1 .. 1000,1 
.. 5); DDASPKproc:-setVars(Y); DDASPKproc:-setInitialTime(T); DDASPKproc:-
setEndTime(TOUT); DDASPKproc:-init(); for i to 999 do retVal := DDASPKproc:-
doIntStep(); Y := DDASPKproc:-getY(); T := DDASPKproc:-getT(); data[i,1] := T;
data[i,2] := Y[1]; data[i,3] := Y[2]; data[i,4] := Y[3]; data[i,5] := Y[4]; if
retVal = 1 then next else break end if end do; imax := i-1; for i to imax do 
if .1e-2 < abs(evalf(data[i,5])-OptValue) then error 
"test 1 using interface restored by RestoreInstance failed"; return false end
if end do; T := 0.; Y := [.5e-1, .6434086, .9702816, 1.250588]; TOUT := 10.0;
data := array(1 .. 1000,1 .. 5); DDASPKproc:-setVars(Y); DDASPKproc:-
setInitialTime(T); DDASPKproc:-setEndTime(TOUT); DDASPKproc:-init(); for i to
999 do retVal := DDASPKproc:-doIntStep(); Y := DDASPKproc:-getY(); T := 
DDASPKproc:-getT(); data[i,1] := T; data[i,2] := Y[1]; data[i,3] := Y[2]; data
[i,4] := Y[3]; data[i,5] := Y[4]; if retVal = 1 then next else break end if 
end do; imax := i-1; if .1e-2 < abs(evalf(data[imax,5])-OptValue) then error 
"test 2 using interface restored by RestoreInstance failed"; return false end
if; Aux:-FileOperations:-removeAllFilesInDir("tmp"); rmdir("tmp"); return true
end proc; simpleProblem1WithAdifor := proc () local g, xRange, yRange, anNLP,
StandardNLP, NPSOLproc, Objf, OptPoint, ExpectedObjf, ExpectedOptPoint, Inform
, aDAESys, StandardModel, DDASPKproc, OptValue, T, Y, TOUT, data, i, retVal, 
imax, result1, result2; ExpectedObjf := 1.00000; ExpectedOptPoint := [1.0, 1.0
]; if Aux:-FileOperations:-dirExists("SimpleProblem1WithAdiforTmp") then error
"temporary directory ./SimpleProblem1WithAdiforTmp needed to run test already \
exists" else mkdir("SimpleProblem1WithAdiforTmp") end if; aDAESys := Aux:-
SystemClasses:-newDAESys(); aDAESys["ODEs"] := [`vm'` = (1-phis)*rhomf/rhom-vm
*qoutoverqin-theta*vm*(lp+ltrm)*Q-2*theta*MWm/rhom*f*kd*zi*If, `vs'` = phis*
rhosf/rhos-vs*qoutoverqin, `zi'` = 1-theta*kd*zi-zi*qoutoverqin, `y'` = rhof*
Cpf/rho/Cp*(yf-yref)-rhof/rho*(y-yref)+negdeltaH/rho/Cp/Tc*theta*vm*lp*rhom/
MWm*Q-theta*alpha*B/rho/Cp*(y-1)]; aDAESys["Parameters"] := [theta = 200.4617,
If = .316531e-1, phis = .6]; aDAESys["DynVars"] := [vm, vs, zi, y]; aDAESys[
"ExplicitAEs"] := [yf = 2/3, yref = 5/9, Tc = 45, T = 45*y, negdeltaH = 21000,
B = 1.4, alpha = 50/3, a3 = -.3495, a2 = -6.7530, a1 = -.4407, f = .8, R = 1.9\
87, MWs = 74.10, MWm = 86.05, phim = 1-phis, rhosf = 777.2814025, rhomf = 918.\
5720, rhomg = 892.0200, rhof = rhomf*phim+rhosf*phis, rhos = 74120/(91.878+.11\
6*T), rhop = 1211-.8496*T, rhom = 958.4-1.3276*T, rhofT = rhom*phim+rhos*phis,
vp = 1-vm-vs, rho = rhom*vm+rhos*vs+rhop*vp, kd = exp(ln(60)+34.99620124-30800
/R/(T+273)), xt = (rhomf*phim-vm*rhomg)/(rhomf*phim+MWm/MWs*phis*rhosf), lt =
exp(a1*xt+a2*xt^2+a3*xt^3), ltrm = 417.3948000*exp(-7569./R/(T+273.)), lp = 
979.8*exp(-4869/R/(T+273))*15^(1/2), Cpf = (.470*rhomf*phim+.716*rhosf*phis)/(
rhomf*phim+rhosf*phis), Cp = (.470*rhom*vm+.716*rhos*vs+rhop*vp*(.321425+.955e\
-3*T))/(rhom*vm+rhos*vs+rhop*vp), Q = 2^(1/2)*(f*kd*zi*If/lt)^(1/2), Rm = (lp+
ltrm)*Q*rhom/MWm*vm+2*f*kd*zi*If, qoutoverqin = rhof/rhofT+theta*MWm*Rm*(1/
rhop-1/rhom)]; aDAESys["AEs"] := []; aDAESys["AlgVars"] := []; Aux:-
SystemClasses:-listOfErrorsInDAESys(aDAESys,'strict'); DDASPKproc := DDASPKSen
:-CreateInstance(aDAESys,"./SimpleProblem1WithAdiforTmp",[theta],'adifor'); 
OptValue := 1.25058770645712380; T := 0.; Y := [.1167175, .6434086, .9702816,
1.250588]; TOUT := 10.0; data := array(1 .. 1000,1 .. 5); DDASPKproc:-setVars(
Y); DDASPKproc:-setInitialTime(T); DDASPKproc:-setEndTime(TOUT); DDASPKproc:-
init(); for i to 999 do retVal := DDASPKproc:-doIntStep(); Y := DDASPKproc:-
getY(); T := DDASPKproc:-getT(); data[i,1] := T; data[i,2] := Y[1]; data[i,3]
:= Y[2]; data[i,4] := Y[3]; data[i,5] := Y[4]; if retVal = 1 then next else 
break end if end do; imax := i-1; result1 := true; for i to imax do if .1e-2 <
abs(evalf(data[i,5])-OptValue) then result1 := false; print(
"Warning(Test1): y[%d]=%f  OptValue=%f",data[i,5],OptValue) end if end do; T 
:= 0.; Y := [.5e-1, .6434086, .9702816, 1.250588]; TOUT := 10.0; data := array
(1 .. 1000,1 .. 5); DDASPKproc:-setVars(Y); DDASPKproc:-setInitialTime(T); 
DDASPKproc:-setEndTime(TOUT); DDASPKproc:-init(); for i to 999 do retVal := 
DDASPKproc:-doIntStep(); Y := DDASPKproc:-getY(); T := DDASPKproc:-getT(); 
data[i,1] := T; data[i,2] := Y[1]; data[i,3] := Y[2]; data[i,4] := Y[3]; data[
i,5] := Y[4]; if retVal = 1 then next else break end if end do; imax := i-1; 
result2 := true; if .1e-2 < abs(evalf(data[imax,5])-OptValue) then result2 :=
false; print("Warning(Test2): y[%d]=%f  OptValue=%f",imx,OptValue) end if; try
dlclose("cwrap_ddaspk1.so"); dlclose("mwrap_cwrapper_ddaspk.dll") catch: NULL
end try; return result1 and result2 end proc; end module; end module;
