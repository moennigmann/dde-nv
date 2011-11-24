NPSOL := module () local Templates, createExternalFunctions, prepareDir, 
createStandardOptionsFile, BigBoundNPSOL, InstanceCounter, OperatingSystem, 
UseAdifor; export CreateInstance, createNumNLP, createNumSysForConstraints, 
createNumSysForCostFunction, createNumSysWrappers, createSharedObject, 
createSharedObjectFromFilelist, getLinConstrOfNLPforNPSOL, 
getConstrBoundsOfNLPforNPSOL, getVarBoundsOfNLPforNPSOL, init, RestoreInstance
, runDefineExternal, SubsIntoSimpleWrapper, Tests; Templates := module () 
export confunTemplate, confunAdiforTemplate, numSysConfunTemplate, 
numSysObjfunTemplate, objfunTemplate, objfunAdiforTemplate, 
SimpleWrapperTemplate; confunTemplate := proc (mode::integer, ncnln::integer,
n::integer, nrowj::integer, needc::integer, x, c, cjac, nstate::integer) local
n_in_confun, ncnln_in_confun, z; global par; declare(par = array(1 .. 
NumParsToBeSubs+1,numeric),x = array(1 .. NumVarsToBeSubs,numeric),c = array(1
.. NumConsToBeSubs,numeric),cjac = array(1 .. NumConsToBeSubs,1 .. 
NumVarsToBeSubs,numeric),n_in_confun = integer,ncnln_in_confun = integer,z = 
array(1 .. NumEAEsToBeSubs+1,numeric)); n_in_confun := NumVarsToBeSubs; if not
n_in_confun = n then ERROR(
`Error in confun: n differs from number of variables`) end if; ncnln_in_confun
:= NumConsToBeSubs; if not ncnln_in_confun = ncnln then ERROR(
`Error in confun: ncnln differs from number of equations`) end if; z := [
EAEsToBeSubs, 0]; c := ConstraintsToBeSubs; cjac := array(1 .. NumConsToBeSubs
,1 .. NumVarsToBeSubs,ListListToBeSubs); RETURN() end proc; 
confunAdiforTemplate := proc (x, c) local z; global par; declare(par = array(1
.. NumParsToBeSubs+1,numeric),x = array(1 .. NumVarsToBeSubs,numeric),c = 
array(1 .. NumConsToBeSubs,numeric),z = array(1 .. NumEAEsToBeSubs+1,numeric))
; z := [EAEsToBeSubs, 0]; c := ConstraintsToBeSubs; RETURN() end proc; 
numSysConfunTemplate := proc (mode::integer, ncnln::integer, n::integer, nrowj
::integer, needc::integer, x, c, cjac, nstate::integer) local i1, i2, g_x, g_c
; global par; declare(par = array(1 .. NumParsToBeSubs+1,numeric),x = array(1
.. NumVarsToBeSubs,numeric),c = array(1 .. NumConsToBeSubs,numeric),cjac = 
array(1 .. NumConsToBeSubs,1 .. NumVarsToBeSubs,numeric),g_x = array(1 .. 
NumVarsToBeSubs,1 .. NumVarsToBeSubs,numeric),g_c = array(1 .. NumVarsToBeSubs
,1 .. NumConsToBeSubs,numeric),i1 = integer,i2 = integer); for i1 to 
NumVarsToBeSubs do for i2 to NumVarsToBeSubs do g_x[i1,i2] := 0. end do; g_x[
i1,i1] := 1.0 end do; DerivsConstraintsAEsNameToBeSubs(n,par,x,g_x,n,c,g_c,n);
ehrpt(); for i1 to NumConsToBeSubs do for i2 to NumVarsToBeSubs do cjac[i1,i2]
:= g_c[i2,i1] end do end do; return end proc; numSysObjfunTemplate := proc (
mode::integer, n::integer, x, objf::float, gradobj, nstate::integer) local i1,
i2, intone, g_x; global par; declare(par = array(1 .. NumParsToBeSubs+1,
numeric),x = array(1 .. NumVarsToBeSubs,numeric),g_x = array(1 .. 
NumVarsToBeSubs,1 .. NumVarsToBeSubs,numeric),gradobj = array(1 .. 
NumVarsToBeSubs,numeric),intone = integer,i1 = integer,i2 = integer); intone 
:= 1; for i1 to NumVarsToBeSubs do for i2 to NumVarsToBeSubs do g_x[i1,i2] :=
0. end do; g_x[i1,i1] := 1.0 end do; DerivsObjectiveAEsNameToBeSubs(n,par,x,
g_x,n,objf,gradobj,n,intone); ehrpt(); return end proc; objfunTemplate := proc
(mode::integer, n::integer, x, objf::float, objgrad, nstate::integer) local 
n_in_objfun, z; global par; declare(par = array(1 .. NumParsToBeSubs+1,numeric
),x = array(1 .. NumVarsToBeSubs,numeric),objgrad = array(1 .. NumVarsToBeSubs
,numeric),n_in_objfun = integer,z = array(1 .. NumEAEsToBeSubs+1,numeric)); 
n_in_objfun := NumVarsToBeSubs; if not n_in_objfun = n then ERROR(
`Error in confun: n differs from number of variables`) end if; z := [
EAEsToBeSubs, 0]; objf := ObjToBeSubs; objgrad := GradToBeSubs; RETURN() end 
proc; objfunAdiforTemplate := proc (x, objf::float) local z; global par; 
declare(par = array(1 .. NumParsToBeSubs+1,numeric),x = array(1 .. 
NumVarsToBeSubs,numeric),z = array(1 .. NumEAEsToBeSubs+1,numeric)); z := [
EAEsToBeSubs, 0]; objf := ObjToBeSubs; RETURN() end proc; end module; 
createExternalFunctions := proc (anNLP::NLP, NameOfConfun::name, NameOfObjfun
::name) local Cons, Pars, Vars, CostFunc, NumConsOfNLP, NumParsOfNLP, 
NumVarsOfNLP, ListListOfJacElements, GradOfCostFunc, confun, objfun, i1, 
ListOfDefExpr, UndefExpr, NumEAEsOfNLP, EAEs, ConsAfterEAEsSubs, 
CostFuncAfterEAEsSubs, pathTemplate; Cons := map(rhs,anNLP["Constraints"]); 
Vars := map(lhs,anNLP["Variables"]); Pars := map(rhs,anNLP["Parameters"]); 
EAEs := map(rhs,anNLP["ExplicitAEs"]); NumConsOfNLP := nops(Cons); 
NumVarsOfNLP := nops(Vars); NumParsOfNLP := nops(Pars); NumEAEsOfNLP := nops(
EAEs); if not Cons = [] and UseAdifor = false then ConsAfterEAEsSubs := Aux:-
ListOperations:-subsEqnListIntoEqn(anNLP["ExplicitAEs"],anNLP["Constraints"]);
ConsAfterEAEsSubs := map(rhs,ConsAfterEAEsSubs); ListListOfJacElements := [seq
([seq(diff(ConsAfterEAEsSubs[i1],Vars[i2]),i2 = 1 .. NumVarsOfNLP)],i1 = 1 ..
NumConsOfNLP)]; NameOfConfun := subs(ParsToBeSubs = op(Pars),
ConstraintsToBeSubs = map(evalf,Cons),EAEsToBeSubs = op(map(evalf,EAEs)),
NumParsToBeSubs = NumParsOfNLP,NumVarsToBeSubs = NumVarsOfNLP,NumConsToBeSubs
= NumConsOfNLP,NumEAEsToBeSubs = NumEAEsOfNLP,ListListToBeSubs = 
ListListOfJacElements,eval(Templates:-confunTemplate)) elif not Cons = [] and
UseAdifor = true then NameOfConfun := subs(ParsToBeSubs = op(Pars),
ConstraintsToBeSubs = map(evalf,Cons),EAEsToBeSubs = op(map(evalf,EAEs)),
NumParsToBeSubs = NumParsOfNLP,NumVarsToBeSubs = NumVarsOfNLP,NumConsToBeSubs
= NumConsOfNLP,NumEAEsToBeSubs = NumEAEsOfNLP,eval(Templates:-
confunAdiforTemplate)) end if; if UseAdifor = false then if not anNLP[
"CostFunction"] = [] then CostFunc := op(anNLP["CostFunction"]); 
GradOfCostFunc := []; CostFuncAfterEAEsSubs := Aux:-ListOperations:-
subsEqnListIntoEqn(anNLP["ExplicitAEs"],CostFunc); for i1 to NumVarsOfNLP do 
GradOfCostFunc := [op(GradOfCostFunc), diff(CostFuncAfterEAEsSubs,Vars[i1])] 
end do; NameOfObjfun := subs(ParsToBeSubs = op(Pars),ObjToBeSubs = evalf(
CostFunc),GradToBeSubs = eval(GradOfCostFunc),EAEsToBeSubs = op(map(evalf,EAEs
)),NumParsToBeSubs = NumParsOfNLP,NumVarsToBeSubs = NumVarsOfNLP,
NumEAEsToBeSubs = NumEAEsOfNLP,eval(Templates:-objfunTemplate)) else 
pathTemplate := anNLP["ObjfunTemplatePath"]; if not Aux:-FileOperations:-
fileExists(pathTemplate) then error "CostFunction in NLP must be not empty or \
ObjfunTemplatePath did not define correct" else read pathTemplate; 
NameOfObjfun := subs(ParsToBeSubs = op(Pars),ObjToBeSubs = evalf(CostFunc),
GradToBeSubs = eval(GradOfCostFunc),EAEsToBeSubs = op(map(evalf,EAEs)),
NumParsToBeSubs = NumParsOfNLP,NumVarsToBeSubs = NumVarsOfNLP,NumEAEsToBeSubs
= NumEAEsOfNLP,eval(objfunTemplate)) end if end if else NameOfObjfun := subs(
ParsToBeSubs = op(Pars),ObjToBeSubs = evalf(CostFunc),EAEsToBeSubs = op(map(
evalf,EAEs)),NumParsToBeSubs = NumParsOfNLP,NumVarsToBeSubs = NumVarsOfNLP,
NumEAEsToBeSubs = NumEAEsOfNLP,eval(Templates:-objfunAdiforTemplate)) end if;
RETURN() end proc; prepareDir := proc (NameOfDir::string) local 
NameOfNPSOLoptionsFile; if not Aux:-FileOperations:-dirExists(NameOfDir) then
mkdir(NameOfDir) end if; NameOfNPSOLoptionsFile := cat(NameOfDir,
"/npsol_options.txt"); if not Aux:-FileOperations:-fileExists(
NameOfNPSOLoptionsFile) then CreateNPSOLstandardOptionsFile(NameOfDir) end if;
if Aux:-FileOperations:-fileExists(cat(NameOfDir,"/mwrap*")) then return false
else return true end if end proc; createStandardOptionsFile := proc (dirname::
string) local fd; fd := fopen(cat("./",dirname,"/npsol_options.txt"),WRITE,
TEXT); fprintf(fd,"begin\n  major iteration limit 9999\n  major print level 20\
\n  Verify level -1\nend\n"); fclose(fd) end proc; BigBoundNPSOL := .1e21; 
InstanceCounter := 0; OperatingSystem := "linux"; CreateInstance := proc (
anNLP::NLP, WorkingDir::string) local npsol, StandardScaledNLP, 
ScaledNLPwithSubstitutions, NpsolInstance, OptArgs, Scaling, ScaledNLP, Ranges
; if 2 < nargs then OptArgs := {args[3 .. -1]} else OptArgs := {} end if; 
UseAdifor := false; if member('adifor',OptArgs) then printf(
"option adifor was requested\n"); UseAdifor := true else printf(
"option adifor was not requested\n") end if; if not prepareDir(WorkingDir) 
then error 
"in working dir, %1, a maple wrapper file named mwrap* already exists", 
WorkingDir end if; if member('scale',OptArgs) then WARNING("SubsIntoSimpleWrap\
per has been changed with respect\nto the modules Scaling it receives"); 
printf("option scale was requested\n"); Scaling := Aux:-NLP:-Scale:-
CreateInstanceForNLP(anNLP); ScaledNLP := Scaling:-subsIntoNLP(anNLP); 
StandardScaledNLP := Aux:-NLP:-subsStandardNotationIntoNLP(ScaledNLP) else 
printf("option scale was not requested\n"); StandardScaledNLP := Aux:-NLP:-
subsStandardNotationIntoNLP(anNLP) end if; if member('scale',OptArgs) then 
StandardScaledNLP := Aux:-NLP:-subsStandardNotationIntoNLP(ScaledNLP) else 
StandardScaledNLP := Aux:-NLP:-subsStandardNotationIntoNLP(anNLP) end if; try
createExternalFunctions(StandardScaledNLP,'npsolconfun','npsolobjfun') catch:
printf("CreateExternalFunctions failed, \n"); printf("  error thrown was %q\n"
,lastexception); printf("  attempting to continue CreateInstance\n") end try;
try npsol := createSharedObject(npsolconfun,npsolobjfun,WorkingDir,
StandardScaledNLP) catch: printf("CreateSharedObject failed, \n"); printf(
"  error thrown was %q\n",lastexception); printf(
"  attempting to continue CreateInstance\n") end try; printf(
"creating instance of interface\n"); if member('scale',OptArgs) then 
ScaledNLPwithSubstitutions := copy(ScaledNLP); ScaledNLPwithSubstitutions[
"Substitutions"] := map(rhs = lhs,StandardScaledNLP["Substitutions"]) else 
ScaledNLPwithSubstitutions := copy(anNLP); ScaledNLPwithSubstitutions[
"Substitutions"] := map(rhs = lhs,StandardScaledNLP["Substitutions"]) end if;
try NpsolInstance := SubsIntoSimpleWrapper(ScaledNLPwithSubstitutions,npsol,
WorkingDir,Scaling,'ScalingPars') catch: printf(
"SubsIntoSimpleWrapper failed, \n"); printf("  error thrown was %q\n",
lastexception); printf("  attempting to continue CreateInstance\n") end try; 
save ScaledNLPwithSubstitutions, cat(WorkingDir,"/NLP.m"); if member('scale',
OptArgs) then Ranges := anNLP["Variables"]; save Ranges, cat(WorkingDir,
"/Ranges.m") end if; RETURN(eval(NpsolInstance)) end proc; createNumNLP := 
proc (Sys::NLP, WorkingDir::string, RestoreKeyword::name, ReqEAEsForObjf::list
(name = term)) local OldDir, ConstraintsAEsInstance, ObjectiveAEsInstance, 
ConstraintsAEsName, ObjectiveAEsName, ConstraintsFileList, Path, 
ObjectiveFileList, FileList, SharedObj, EAEsForObjf; if not member(
RestoreKeyword,{'new', 'restore'}) then error 
"3rd input argument, RestoreKeyword, must be new or restore" end if; if not 
type(NumSys,`module`) then printf("reading module NumSys\n"); read cat(
_EnvModulesDir,"/NumSys/source/NumSys.mpl") end if; if not type(ADIFOR,
`module`) then printf("reading module ADIFOR\n"); read cat(_EnvModulesDir,
"/ADIFOR/source/ADIFOR.mpl") end if; if not Aux:-FileOperations:-dirExists(
WorkingDir) then system(cat("mkdir ",WorkingDir)) end if; OldDir := currentdir
(); if 3 < nargs then EAEsForObjf := ReqEAEsForObjf else EAEsForObjf := '
EAEsForObjf' end if; try ConstraintsAEsInstance := NPSOL:-
createNumSysForConstraints(Sys,WorkingDir,RestoreKeyword) catch: printf(
"creation of EAEs and AEs for constraints failed with error:\n"); error 
finally currentdir(OldDir) end try; try currentdir(WorkingDir); system(
"mkdir CreateNumSysForCostFunction"); ObjectiveAEsInstance := NPSOL:-
createNumSysForCostFunction(Sys,"CreateNumSysForCostFunction",RestoreKeyword,
EAEsForObjf) catch: printf(
"creation of EAEs and AEs for objective function failed with error:\n"); error
finally currentdir(OldDir) end try; ConstraintsAEsName := convert(cat("g_",
ConstraintsAEsInstance:-GetFileList()[1]),name); ObjectiveAEsName := convert(
cat("g_",ObjectiveAEsInstance:-GetFileList()[1]),name); try currentdir(
WorkingDir); NPSOL:-createNumSysWrappers(nops(Sys["Parameters"]),nops(Sys[
"Variables"]),nops(Sys["Constraints"]),ConstraintsAEsName,ObjectiveAEsName) 
catch: printf("error generating wrappers:\n"); error finally currentdir(OldDir
) end try; ConstraintsFileList := ConstraintsAEsInstance:-GetFileList(); Path
:= ConstraintsAEsInstance:-GetHomeDir(); ConstraintsFileList := [seq(cat(Path,
"/g_",ConstraintsFileList[i1]),i1 = 1 .. nops(ConstraintsFileList))]; 
ObjectiveFileList := ObjectiveAEsInstance:-GetFileList(); Path := 
ObjectiveAEsInstance:-GetHomeDir(); ObjectiveFileList := [seq(cat(Path,"/g_",
ObjectiveFileList[i1]),i1 = 1 .. nops(ObjectiveFileList))]; FileList := [
"npsolconfun", "npsolobjfun", op(ConstraintsFileList), op(ObjectiveFileList)];
SharedObj := NPSOL:-createSharedObjectFromFilelist(FileList,WorkingDir); 
return SharedObj end proc; createNumSysForConstraints := proc (Sys::NLP, 
WorkingDir::string, RestoreKeyword) local OldDir, EAEsInstance, VisEAEs, 
ConstraintsAEsInstance; if not member(RestoreKeyword,{'new', 'restore'}) then
error "3rd input argument, RestoreKeyword, must be new or restore" end if; if
not type(NumSys,`module`) then printf("reading module NumSys\n"); read cat(
_EnvModulesDir,"/NumSys/source/NumSys.mpl") end if; if not type(ADIFOR,
`module`) then printf("reading module ADIFOR\n"); read cat(_EnvModulesDir,
"/ADIFOR/source/ADIFOR.mpl") end if; if not Aux:-FileOperations:-dirExists(
WorkingDir) then system(cat("mkdir ",WorkingDir)) end if; OldDir := currentdir
(); try EAEsInstance := NumSys:-EAE:-CreateInstance(Sys,WorkingDir,
RestoreKeyword,'NoInterface') catch: printf(
"error while creating instance of NumSys for EAEs:\n"); error finally 
currentdir(OldDir) end try; VisEAEs := map(lhs,Sys[ExplicitAEs]); try 
ConstraintsAEsInstance := NumSys:-AE:-CreateInstance(Sys,WorkingDir,
RestoreKeyword,'NoInterface',cat(WorkingDir,
"/NumSys-EAE-CreateInstanceForStandardEAEs"),VisEAEs) catch: printf("in Create\
NumSysForConstraints, error while calling NumSys:-AE:-CreateInstance:\n"); 
error finally currentdir(OldDir) end try; try currentdir(WorkingDir); ADIFOR:-
RunAD(ConstraintsAEsInstance:-GetFileList(),'y',nops(Sys[Variables]),'res') 
catch: printf(
"in CreateNumSysForConstraints, error while calling ADIFOR:-RunAD:\n"); error
finally currentdir(OldDir) end try; return ConstraintsAEsInstance end proc; 
createNumSysForCostFunction := proc (Sys::NLP, WorkingDir::string, 
RestoreKeyword::name, EAEsForObjf::{name, list(name = term), []}) local 
CostFuncNLP, ObjectiveAEsInstance, OldDir; OldDir := currentdir(); if Sys[
"CostFunction"] = [] then error "cost function is empty" end if; CostFuncNLP 
:= copy(Sys); CostFuncNLP["Constraints"] := [0 = op(CostFuncNLP["CostFunction"
])]; CostFuncNLP["CostFunction"] := [0]; if not type(EAEsForObjf,name) then 
CostFuncNLP["ExplicitAEs"] := EAEsForObjf end if; try ObjectiveAEsInstance :=
createNumSysForConstraints(CostFuncNLP,WorkingDir,RestoreKeyword) catch: 
printf("in CreateNumSysForCostFunction, call to CreateNumSysForConstraints fai\
led with following error:\n"); error finally currentdir(OldDir) end try; 
return ObjectiveAEsInstance end proc; createNumSysWrappers := proc (
NumParsInNLP::{0, posint}, NumVarsInNLP::posint, NumConsInNLP::posint, 
DerivsConstraintsAEsName::name, DerivsObjectiveAEsName::name) local 
npsolconfun, npsolobjfun; npsolconfun := subs(NumParsToBeSubs = NumParsInNLP,
NumVarsToBeSubs = NumVarsInNLP,NumConsToBeSubs = NumConsInNLP,
DerivsConstraintsAEsNameToBeSubs = DerivsConstraintsAEsName,eval(Templates:-
numSysConfunTemplate)); system("rm -f npsolconfun.f"); codegen[fortran](
npsolconfun,filename = "npsolconfun.f",mode = double); npsolobjfun := subs(
NumParsToBeSubs = NumParsInNLP,NumVarsToBeSubs = NumVarsInNLP,NumConsToBeSubs
= NumConsInNLP,DerivsObjectiveAEsNameToBeSubs = DerivsObjectiveAEsName,eval(
Templates:-numSysObjfunTemplate)); system("rm -f npsolobjfun.f"); codegen[
fortran](npsolobjfun,filename = "npsolobjfun.f",mode = double); return end 
proc; createSharedObject := proc (NpsolConfun::name, NpsolObjfun::name, 
NpsolDir::string, StandardNLP) local oldDir, LinkerCommand, npsol, 
systemsuccessful, N, NCLIN, NCNLN, NROWA, NROWJ, NROWR, A, BL, BU, INFORM, 
ITER, ISTATE, C, CJAC, CLAMBDA, OBJF, GRAD, R, XVEC, LENIW, LENW, 
systemsuccessful2, SharedObjectName, fd1, fd2; InstanceCounter := 
InstanceCounter+1; oldDir := currentdir(); currentdir(NpsolDir); if Aux:-
FileOperations:-fileExists("npsolconfun.f") then currentdir(oldDir); error 
"routine npsolconfun.f already exists in directory %1", NpsolDir end if; if 
Aux:-FileOperations:-fileExists("npsolobjfun.f") then currentdir(oldDir); 
error "routine npsolobjfun.f already exists in directory %1", NpsolDir end if;
if not type(eval(NpsolConfun),name) then codegen[fortran](NpsolConfun,filename
= "npsolconfun.f",precision = double); if UseAdifor = true then ADIFOR:-RunAD(
"npsolconfun",['x'],nops(StandardNLP["Variables"]),'c'); fd1 := fopen(cat(
"./sedcommandforconfun.txt"),WRITE,TEXT); fprintf(fd1,
"sed -e 's/g_npsolconfun/gnpsolconfun/g' g_npsolconfun.f > gnpsolconfun.f"); 
fclose(fd1); system("sh sedcommandforconfun.txt"); system(
"rm -f sedcommandforconfun.txt"); system("rm -f g_npsolconfun.f") end if; if 
OperatingSystem = "linux" and UseAdifor = false then system(
"gfortran -c -fPIC npsolconfun.f -o npsolconfun.o") elif OperatingSystem = 
"linux" and UseAdifor = true then system(
"gfortran -c -fPIC gnpsolconfun.f -o gnpsolconfun.o") end if else if 
OperatingSystem = "linux" then systemsuccessful := system(cat("cp ",
_ModulesDirectory,"/Npsol/ext_routines/f2c/dummies/npsolconfun.o .")); 
systemsuccessful2 := system(cat("cp ",_ModulesDirectory,
"/Npsol/ext_routines/f2c/dummies/gnpsolconfun.o .")) end if; if not 
systemsuccessful = 0 then error "dummy file npsolconfun.o not found" end if; 
if not systemsuccessful2 = 0 then error "dummy file gnpsolconfun.o not found"
end if end if; codegen[fortran](NpsolObjfun,filename = "npsolobjfun.f",
precision = double); if UseAdifor = true then ADIFOR:-RunAD("npsolobjfun",['x'
],nops(StandardNLP["Variables"]),'objf'); fd2 := fopen(cat(
"./sedcommandforobjfun.txt"),WRITE,TEXT); fprintf(fd2,
"sed -e 's/g_npsolobjfun/gnpsolobjfun/g' g_npsolobjfun.f > gnpsolobjfun.f"); 
fclose(fd2); system("sh sedcommandforobjfun.txt"); system(
"rm -f sedcommandforobjfun.txt"); system("rm -f g_npsolobjfun.f") end if; if 
OperatingSystem = "linux" and UseAdifor = false then system(
"gfortran -c -fPIC npsolobjfun.f -o npsolobjfun.o") elif OperatingSystem = 
"linux" and UseAdifor = true then system(
"gfortran -c -fPIC gnpsolobjfun.f -o gnpsolobjfun.o") end if; SharedObjectName
:= cat("cwrap_npsol",InstanceCounter,".so"); if OperatingSystem = "linux" and
UseAdifor = false then LinkerCommand := cat(
`gcc -shared -Xlinker -Bsymbolic  npsolconfun.o npsolobjfun.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/npsol_open.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/npsol_close.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/set_params.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/cwrap_npsol.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/shared_obj/npsol.so `,
`-lf2c -lm -lgfortran -o `,SharedObjectName) elif OperatingSystem = "linux" 
and UseAdifor = true then LinkerCommand := cat(
`gcc -shared -Xlinker -Bsymbolic gnpsolconfun.o gnpsolobjfun.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/npsolconfunad.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/npsolobjfunad.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/npsol_open.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/npsol_close.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/set_params.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/cwrap/cwrap_npsol.o `,
_ModulesDirectory,`/Npsol/ext_routines/f2c/shared_obj/npsol.so `,_AD_LIB,
`/lib/ReqADIntrinsics-Linux86_64.o `,_AD_LIB,
`/lib/libADIntrinsics-Linux86_64.a `,`-lf2c -lm -lgfortran -o `,
SharedObjectName) end if; system(LinkerCommand); npsol := runDefineExternal(
SharedObjectName); currentdir(oldDir); RETURN(npsol) end proc; 
createSharedObjectFromFilelist := proc (FileList::list({name, string}), 
WorkingDir::string) local SharedObjectName, LinkerCommand, OldDir, 
ObjectListForLinker, NpsolCwrapper, StringFileList, FortranFileList, 
ObjectList, item, CompileCommand, Path, File; InstanceCounter := 
InstanceCounter+1; OldDir := currentdir(); SharedObjectName := cat(
"cwrap_npsol",InstanceCounter,".so"); StringFileList := map(convert,FileList,
string); ObjectList, FortranFileList := Aux:-FileOperations:-findObjectFile(
StringFileList); FortranFileList := [seq(cat(FortranFileList[i1],".f"),i1 = 1
.. nops(FortranFileList))]; if OperatingSystem = "linux" then try currentdir(
WorkingDir); for item in FortranFileList while true do Path, File := Aux:-
FileOperations:-splitPathAndFile(item); Aux:-FileOperations:-
runSystemCommandInDir(cat("g77 -O0 -c ",File),Path); ObjectList := [op(
ObjectList), cat(Path,File[1 .. -3],".o")] end do; ObjectListForLinker := cat(
seq(cat(ObjectList[i1]," "),i1 = 1 .. nops(ObjectList))) catch: printf(
"building of file lists failed with error:\n"); error finally currentdir(
OldDir) end try; LinkerCommand := cat("gcc -O0 -shared -Xlinker -Bsymbolic ",
ObjectListForLinker,_EnvModulesDir,
`/NPSOL/ext_routines/f2c/cwrap/npsol_open.o `,_EnvModulesDir,
`/NPSOL/ext_routines/f2c/cwrap/npsol_close.o `,_EnvModulesDir,
`/NPSOL/ext_routines/f2c/cwrap/set_params.o `,_EnvModulesDir,
"/NPSOL/ext_routines/f2c/cwrap/cwrap_npsol.o ",_EnvModulesDir,
"/NPSOL/ext_routines/f2c/shared_obj/npsol.so ",
"$AD_LIB/lib/ReqADIntrinsics-$AD_OS.o ",
"$AD_LIB/lib/libADIntrinsics-$AD_OS.a ","$AD_LIB/lib/libSparsLinC-$AD_OS.a ",
"-lf2c -lm -o ",SharedObjectName) else error "operating system unknown" end if
; try currentdir(WorkingDir); printf("%s\n",LinkerCommand); system(
LinkerCommand) catch: printf("linking failed with following error: "); error 
finally currentdir(OldDir) end try; try currentdir(WorkingDir); NpsolCwrapper
:= runDefineExternal(SharedObjectName) catch: printf(
"define_external failed with following error"); error finally currentdir(
OldDir) end try; return NpsolCwrapper end proc; getLinConstrOfNLPforNPSOL := 
proc (anNLP::NLP, NameForAMAT::name, NameForbVEC::name) local AMAT, bVEC, 
NCLIN, LinConstrParsAndEAEsSubstituted, N, LinConstrEAEsSubstituted; if not 
anNLP["LinearConstraints"] = [] then NCLIN := nops(anNLP["LinearConstraints"])
; LinConstrEAEsSubstituted := Aux:-ListOperations:-subsEqnListIntoEqn(anNLP[
"ExplicitAEs"],anNLP["LinearConstraints"]); LinConstrParsAndEAEsSubstituted :=
subs(anNLP["Parameters"],map(rhs,LinConstrEAEsSubstituted)); Aux:-LinearEqns:-
checkLinearityOfIn(LinConstrParsAndEAEsSubstituted,map(lhs,anNLP["Variables"])
,'AMAT','bVEC') else NCLIN := 0; N := nops(anNLP["Variables"]); AMAT := matrix
(1,N,[seq(0,i1 = 1 .. N)]) end if; NameForAMAT := AMAT; NameForbVEC := bVEC; 
return end proc; getConstrBoundsOfNLPforNPSOL := proc (anNLP::NLP, bVEC::{name
, vector}, NameForLowerBounds::name, NameForUpperBounds::name) local NCLIN, 
NCNLN, LowerBoundsConstraints, UpperBoundsConstraints, i1, item; NCLIN := nops
(anNLP["LinearConstraints"]); NCNLN := nops(anNLP["Constraints"]); 
LowerBoundsConstraints := [seq(bVEC[i1],i1 = 1 .. NCLIN), seq(0,i1 = NCLIN+1 
.. NCLIN+NCNLN)]; UpperBoundsConstraints := []; if not anNLP[
"LinearConstraints"] = [] then for i1 to nops(anNLP["LinearConstraints"]) do 
item := anNLP["LinearConstraints"][i1]; if type(item,`=`) then 
UpperBoundsConstraints := [op(UpperBoundsConstraints), bVEC[i1]] elif type(
item,0 < term) then UpperBoundsConstraints := [op(UpperBoundsConstraints), 
BigBoundNPSOL] else error 
"internal error, this only works if inequality has standard form 0 < expr." 
end if end do end if; for i1 to nops(anNLP["Constraints"]) do item := anNLP[
"Constraints"][i1]; if type(item,`=`) then UpperBoundsConstraints := [op(
UpperBoundsConstraints), 0] elif type(item,0 < term) then 
UpperBoundsConstraints := [op(UpperBoundsConstraints), BigBoundNPSOL] else 
error "internal error, nonlinear constraints are expected to have form 0=... o\
r 0<..." end if end do; NameForUpperBounds := UpperBoundsConstraints; 
NameForLowerBounds := LowerBoundsConstraints; return end proc; 
getVarBoundsOfNLPforNPSOL := proc (anNLP::NLP, NameForLowerBounds::name, 
NameForUpperBounds::name) local RangesVars, LowerBoundsVars, UpperBoundsVars;
if ScalingOn = true then RangesVars := Aux:-NLP:-Scale:-
replaceInfiniteRangesByTrivialRanges(anNLP["Variables"]) end if; RangesVars :=
map(rhs,anNLP["Variables"]); NameForLowerBounds := [seq(subs(infinity = 
BigBoundNPSOL,-infinity = -BigBoundNPSOL,op(1,RangesVars[i1])),i1 = 1 .. nops(
RangesVars))]; NameForUpperBounds := [seq(subs(infinity = BigBoundNPSOL,-
infinity = -BigBoundNPSOL,op(2,RangesVars[i1])),i1 = 1 .. nops(RangesVars))];
return end proc; init := proc () local OperatingSystemTmp, 
OperatingSystemTmpSol, OperatingSystemTmpWin, OperatingSystemTmpLin; 
InstanceCounter := 0; OperatingSystemTmpSol := getenv(OSTYPE); 
OperatingSystemTmpWin := getenv(OS); OperatingSystemTmpLin := getenv(OSNAME);
OperatingSystemTmp := cat(OperatingSystemTmpSol,OperatingSystemTmpWin,
OperatingSystemTmpLin); if not searchtext("solaris",OperatingSystemTmp) = 0 
then OperatingSystemTmp := "solaris" end if; if not searchtext("windows",
OperatingSystemTmp) = 0 then OperatingSystemTmp := "windows" end if; if not 
searchtext("inux",OperatingSystemTmp) = 0 then OperatingSystemTmp := "linux" 
end if; if not `in`(OperatingSystemTmp,{"linux", "solaris", "windows"}) then 
error "OperatingSystem could not be detected. Set environment variable OSTYPE,\
 OS, or OSNAME!" end if; OperatingSystem := OperatingSystemTmp end proc; 
RestoreInstance := proc (WorkingDir::string) local NPSOLinstance, 
SharedObjectName, FileName, SimpleWrapper, OldDir, NumPars, StandardSys, npsol
, Scaling, VariablesTmp; OldDir := currentdir(); currentdir(WorkingDir); 
FileName := "NLP.m"; if not Aux:-FileOperations:-fileExists(FileName) then 
currentdir(OldDir); error 
"in the requested dir %1, file named NLP.m does not exist", WorkingDir end if;
read "NLP.m"; if 1 < nargs then if not type(args[2],string) then error 
"optional third argument, SharedObjectName, must be a string" else 
SharedObjectName := args[2] end if else SharedObjectName := Aux:-
FileOperations:-findFirstSharedObject("cwrap_npsol") end if; if not Aux:-
FileOperations:-fileExists(SharedObjectName) then currentdir(OldDir); error 
"in the requested dir %1, file named %2 does not exist", WorkingDir, 
SharedObjectName end if; if not 1 < nargs then warning(
"using shared object %1 in %2",SharedObjectName,WorkingDir) end if; npsol := 
define_external('cwrapper_npsol',N::integer[4],NCLIN::integer[4],NCNLN::
integer[4],NROWA::integer[4],NROWJ::integer[4],NROWR::integer[4],A::ARRAY(1 ..
NROWA,1 .. N,float[8]),B::ARRAY(1 .. 3,1 .. N+NCLIN+NCNLN,float[8]),INFORM::
REF(integer[4]),ITER::REF(integer[4]),ISTATE::ARRAY(1 .. N+NCLIN+NCNLN,integer
[4]),C::ARRAY(1 .. NCNLN,float[8]),CJAC::ARRAY(1 .. NROWJ,1 .. N,float[8]),
OBJF::REF(float[8]),GRAD::REF(ARRAY(1 .. N,float[8])),R::ARRAY(1 .. NROWR,1 ..
N,float[8]),XVEC::ARRAY(1 .. N,float[8]),LE::ARRAY(1 .. 2,integer[4]),
NumParams::integer[4],Params::ARRAY(1 .. NumParams,float[8]),LIB = 
SharedObjectName); if Aux:-FileOperations:-fileExists("Ranges.m") then read 
"Ranges.m"; Scaling := Aux:-NLP:-Scale:-CreateInstanceForNLP(Ranges) end if; 
currentdir(OldDir); NPSOLinstance := SubsIntoSimpleWrapper(
ScaledNLPwithSubstitutions,npsol,WorkingDir,Scaling,'ScalingPars'); currentdir
(OldDir); unassign('ScaledNLPwithSubstitutions'); unassign('Ranges'); RETURN(
eval(NPSOLinstance)) end proc; runDefineExternal := proc (SharedObjectName::
string) local npsol; npsol := define_external('cwrapper_npsol',N::integer[4],
NCLIN::integer[4],NCNLN::integer[4],NROWA::integer[4],NROWJ::integer[4],NROWR
::integer[4],A::ARRAY(1 .. NROWA,1 .. N,float[8]),B::ARRAY(1 .. 3,1 .. N+NCLIN
+NCNLN,float[8]),INFORM::REF(integer[4]),ITER::REF(integer[4]),ISTATE::ARRAY(1
.. N+NCLIN+NCNLN,integer[4]),C::ARRAY(1 .. NCNLN,float[8]),CJAC::ARRAY(1 .. 
NROWJ,1 .. N,float[8]),OBJF::REF(float[8]),GRAD::ARRAY(1 .. N,float[8]),R::
ARRAY(1 .. NROWR,1 .. N,float[8]),XVEC::ARRAY(1 .. N,float[8]),LE::ARRAY(1 ..
2,integer[4]),NumParams::integer[4],Params::ARRAY(1 .. NumParams,float[8]),LIB
= SharedObjectName); return npsol end proc; SubsIntoSimpleWrapper := proc (
anNLP::NLP, NpsolCwrapper::procedure, WorkingDir::string, ReqScaling, 
ReqScalingPars) local LowerBoundsVars, UpperBoundsVars, LowerBoundsConstr, 
UpperBoundsConstr, Amat, bVEC, NCNLN, NROWA, NROWJ, NROWR, LENIW, N, NCLIN, 
LENW, SimpleWrapperInstance, BL_swnpsol, BU_swnpsol, xvec_swnpsol, A_swnpsol,
istate_swnpsol, c_swnpsol, cjac_swnpsol, clambda_swnpsol, grad_swnpsol, 
r_swnpsol, SimpleWrapperTemplate, RequiredVars, ParNames, NpsolOptions, 
NumParams, Params, KnownNpsolOptions; N := nops(anNLP["Variables"]); NCLIN :=
nops(anNLP["LinearConstraints"]); getLinConstrOfNLPforNPSOL(anNLP,Amat,bVEC);
NCNLN := nops(anNLP["Constraints"]); NROWA := max(1,NCLIN); NROWJ := max(1,
NCNLN); NROWR := N; LENIW := 3*N+NCLIN+2*NCNLN; if NCLIN = 0 and NCNLN = 0 
then LENW := 20*N elif NCNLN = 0 and NCNLN <> 0 then LENW := 2*N^2+20*N+11*
NCLIN else LENW := 2*N^2+N*NCLIN+2*N*NCNLN+20*N+11*NCLIN+21*NCNLN end if; 
getVarBoundsOfNLPforNPSOL(anNLP,LowerBoundsVars,UpperBoundsVars); 
getConstrBoundsOfNLPforNPSOL(anNLP,bVEC,LowerBoundsConstr,UpperBoundsConstr);
BL_swnpsol := Vector(1 .. N+NCLIN+NCNLN,map(evalf,[op(LowerBoundsVars), op(
LowerBoundsConstr)]),datatype = float[8],order = C_order); BU_swnpsol := 
Vector(1 .. N+NCLIN+NCNLN,map(evalf,[op(UpperBoundsVars), op(UpperBoundsConstr
)]),datatype = float[8],order = C_order); A_swnpsol := Matrix(1 .. NROWA,1 ..
N,map(evalf,Amat),datatype = float[8],order = C_order); NumParams := nops(
anNLP["Parameters"]); Params := Vector(1 .. NumParams,map(rhs,anNLP[
"Parameters"]),datatype = float[8]); xvec_swnpsol := Vector(1 .. N,datatype =
float[8],order = C_order); istate_swnpsol := Vector(1 .. N+NCLIN+NCNLN,
datatype = integer[4]); c_swnpsol := Vector(1 .. NCNLN,datatype = float[8]); 
cjac_swnpsol := Matrix(1 .. NROWJ,1 .. N,datatype = float[8],order = C_order);
clambda_swnpsol := Vector(1 .. N+NCLIN+NCNLN,datatype = float[8]); 
grad_swnpsol := Vector(1 .. N,datatype = float[8]); r_swnpsol := Matrix(1 .. 
NROWR,1 .. N,datatype = float[8],order = C_order); KnownNpsolOptions := {
"verify level", "derivative level", "function precision", "linear feasibility"
, "major iteration limit", "major print level", "minor iteration limit", 
"minor print level", "nonlinear feasibility", "optimality tolerance"}; 
NpsolOptions := table(); NpsolOptions["major iteration limit"] := 999; 
NpsolOptions["major print level"] := 20; if anNLP["CostFunction"] = [] then 
NpsolOptions["derivative level"] := 2 end if; RequiredVars := map(lhs,anNLP[
"Variables"]); ParNames := map(lhs,anNLP["Parameters"]); SimpleWrapperTemplate
:= module () local TheNLP, Scaling, ScalingOn, ScalingPars, ScalingParsOn; 
export boundsNotViolated, createDebugFiles, createOptionsFile, 
getLowerActiveVarBounds, getUpperActiveVarBounds, getBL, getBU, getC, getCjac,
getClambda, getObjf, getGrad, getIndicesActiveVarBounds, 
getIndicesLinearInequalityConstraints, 
getIndicesNonlinearInequalityConstraints, getInform, getIter, getIstate, getR,
getOptions, getOptPoint, getParameters, getSys, getVarBounds, getVariables, 
getXVEC, runOpt, setAmat, setBL, setBU, setIstate, setCjac, setClambda, setR,
setLowerBoundsConstr, setUpperBoundsConstr, setLowerBoundsForVar, 
setUpperBoundsForVar, setLowerBoundsVars, setUpperBoundsVars, setBoundsForVar,
setOption, setParameters, setVariables, setXVEC; TheNLP := copy(anNLP); if 
type(ReqScaling,`module`) then Scaling := copy(ReqScaling); ScalingOn := true
else ScalingOn := false end if; if type(ReqScalingPars,`module`) then 
ScalingPars := copy(ReqScalingPars); ScalingParsOn := true else ScalingParsOn
:= false end if; boundsNotViolated := proc () local ViolatedLowerBounds, 
ViolatedUpperBounds, i1, NamesViolated; ViolatedLowerBounds := []; 
ViolatedUpperBounds := []; for i1 to N do if xvec_swnpsol[i1] < BL_swnpsol[i1]
then ViolatedLowerBounds := [op(ViolatedLowerBounds), i1] end if; if 
BU_swnpsol[i1] < xvec_swnpsol[i1] then ViolatedUpperBounds := [op(
ViolatedUpperBounds), i1] end if end do; if not ViolatedUpperBounds = [] then
NamesViolated := [seq(RequiredVars[i1],i1 = ViolatedUpperBounds)]; WARNING(
"some variables violate upper bounds: %1",NamesViolated) end if; if not 
ViolatedLowerBounds = [] then NamesViolated := [seq(RequiredVars[i1],i1 = 
ViolatedLowerBounds)]; WARNING("some variables violate lower bounds: %1",
NamesViolated) end if; if not (ViolatedLowerBounds = [] and 
ViolatedUpperBounds = []) then return false else return true end if end proc;
createDebugFiles := proc () local fd, i1, dirname; if 0 < nargs then if type(
arg[1] = string) then dirname := arg[1] else error 
"expecting first input argument to be of type string" end if else dirname := 
cat(_ModulesDirectory,"/Npsol/ext_routines/debugging/") end if; printf(
"Writing debug files to directory %s",dirname); fd := fopen(cat(dirname,
"/system_dimension.h"),WRITE,TEXT); fprintf(fd,"ncnln= %d;\n",NCNLN); fprintf(
fd,"n    = %d;\n",N); fprintf(fd,"nrowj= %d;\n",NROWJ); fprintf(fd,
"numparams= %d;\n",NumParams); fclose(fd); fd := fopen(cat(dirname,
"/variables.h"),WRITE,TEXT); for i1 to N do fprintf(fd,
"x[%d]= (double) %e; /* %s */\n",i1-1,xvec_swnpsol[i1],convert(RequiredVars[i1
],string)) end do; fprintf(fd,"\n"); for i1 to NumParams do fprintf(fd,
"par[%d]= (double) %e; // %s\n",i1-1,Params[i1],convert(ParNames[i1],string))
end do; fclose(fd) end proc; createOptionsFile := proc () local item, 
EntryNames, OldDir, fd; OldDir := currentdir(); currentdir(WorkingDir); fd :=
fopen(cat("./npsol_options.txt"),WRITE,TEXT); fprintf(fd,"begin\n"); 
EntryNames := map(op,[indices(NpsolOptions)]); for item in EntryNames while 
true do fprintf(fd,cat("  ",item," ",convert(NpsolOptions[item],string),"\n"))
end do; fprintf(fd,"end\n"); fclose(fd); currentdir(OldDir); return end proc;
getBL := proc () return eval(BL_swnpsol) end proc; getBU := proc () return 
eval(BU_swnpsol) end proc; getInform := proc () return inform_swnpsol end proc
; getIter := proc () return iter_swnpsol end proc; getIstate := proc () return
istate_swnpsol end proc; getC := proc () return c_swnpsol end proc; getCjac :=
proc () return cjac_swnpsol end proc; getClambda := proc () return 
clambda_swnpsol end proc; getObjf := proc () return objf_swnpsol end proc; 
getGrad := proc () return grad_swnpsol end proc; getR := proc () return 
r_swnpsol end proc; getOptions := proc () local KnownOptions, item; 
KnownOptions := map(op,[indices(NpsolOptions)]); printf(
"current options are:\n"); for item in KnownOptions while true do printf(
"%s= %d\n",convert(item,string),NpsolOptions[item]) end do; printf(
"known options are:\n"); for item in KnownNpsolOptions while true do printf(
"%s, ",convert(item,string)) end do; printf("\n"); return end proc; 
getParameters := proc () local current, RescaledCurrent, c, InvD, i1; if 
ScalingParsOn = true then RescaledCurrent := array(1 .. NumParams); InvD := 
ScalingPars:-getInvD(); c := ScalingPars:-getc(); for i1 to NumParams do 
RescaledCurrent[i1] := ParNames[i1] = InvD[i1]*(Params[i1]+c[i1]) end do; 
RescaledCurrent := convert(RescaledCurrent,list); return RescaledCurrent else
current := [seq(ParNames[i1] = Params[i1],i1 = 1 .. NumParams)]; return 
current end if end proc; getSys := proc () WARNING(
"scaled system is returned if scaling was requested"); return eval(anNLP) end
proc; getVariables := proc () local NamedList, OptPoint; OptPoint := getXVEC()
; NamedList := [seq(RequiredVars[i1] = OptPoint[i1],i1 = 1 .. nops(
RequiredVars))]; return NamedList end proc; getXVEC := proc () local InvD, c,
i1, RescaledXVEC; if ScalingOn = true then RescaledXVEC := array(1 .. N); InvD
:= Scaling:-getInvD(); c := Scaling:-getc(); for i1 to N do RescaledXVEC[i1] 
:= InvD[i1]*(xvec_swnpsol[i1]+c[i1]) end do; return eval(RescaledXVEC) else 
return eval(xvec_swnpsol) end if end proc; setIstate := proc (ReqVal::list(
integer)) istate_swnpsol := Vector(1 .. N+NCLIN+NCNLN,ReqVal,datatype = 
integer[4],order = C_order); return end proc; setCjac := proc (ReqVal::matrix)
cjac_swnpsol := Matrix(1 .. NROWJ,1 .. N,map(evalf,ReqVal),datatype = float[8]
,order = Fortran_order); return end proc; setClambda := proc (ReqVal::list(
EvalsToFloat)) clambda_swnpsol := Vector(1 .. N+NCLIN+NCNLN,map(evalf,ReqVal),
datatype = float[8],order = C_order); return end proc; setR := proc (ReqVal::
matrix) r_swnpsol := Matrix(1 .. NROWR,1 .. N,map(evalf,ReqVal),datatype = 
float[8],order = Fortran_order); return end proc; setBL := proc () BL_swnpsol
:= Vector(1 .. N+NCLIN+NCNLN,map(evalf,[op(LowerBoundsVars), op(
LowerBoundsConstr)]),datatype = float[8],order = C_order); return end proc; 
setBU := proc () BU_swnpsol := Vector(1 .. N+NCLIN+NCNLN,map(evalf,[op(
UpperBoundsVars), op(UpperBoundsConstr)]),datatype = float[8],order = C_order)
; return end proc; setAmat := proc (ReqVal::matrix) A_swnpsol := Matrix(1 .. 
NROWA,1 .. N,map(evalf,ReqVal),datatype = float[8],order = Fortran_order); 
return end proc; setLowerBoundsConstr := proc (ReqVal::list(EvalsToFloat)) if
not nops(ReqVal) = NCLIN+NCNLN then error 
"expecting %1 values for lower bounds of variables", N end if; 
LowerBoundsConstr := ReqVal; setBL(); return end proc; setUpperBoundsConstr :=
proc (ReqVal::list(EvalsToFloat)) if not nops(ReqVal) = NCLIN+NCNLN then error
"expecting %1 values for lower bounds of variables", N end if; 
UpperBoundsConstr := ReqVal; setBU(); return end proc; setLowerBoundsVars := 
proc (ReqVal::list(EvalsToFloat)) if not nops(ReqVal) = N then error 
"expecting %1 values for lower bounds of variables", N end if; LowerBoundsVars
:= ReqVal; setBL(); return end proc; setUpperBoundsVars := proc (ReqVal::list(
EvalsToFloat)) if not nops(ReqVal) = N then error 
"expecting %1 values for lower bounds of variables", N end if; UpperBoundsVars
:= ReqVal; setBU(); return end proc; runOpt := proc () local NextNumber, 
NextFileName, OldDir, B_swnpsol, LE_swnpsol; OldDir := currentdir(); 
currentdir(WorkingDir); B_swnpsol := Matrix([[convert(BL_swnpsol,array)], [
convert(BU_swnpsol,array)], [convert(clambda_swnpsol,array)]],datatype = float
[8],order = C_order); LE_swnpsol := Vector([LENIW, LENW],datatype = integer[4]
); NpsolCwrapper(N,NCLIN,NCNLN,NROWA,NROWJ,NROWR,A_swnpsol,B_swnpsol,'
inform_swnpsol','iter_swnpsol',istate_swnpsol,c_swnpsol,cjac_swnpsol,'
objf_swnpsol',grad_swnpsol,r_swnpsol,xvec_swnpsol,LE_swnpsol,NumParams,Params)
; NextNumber := Aux:-FileOperations:-nextFileNumber("NpsolOut"); NextFileName
:= cat("NpsolOut.",NextNumber); Aux:-FileOperations:-renameFile("fort.9",
NextFileName); currentdir(OldDir); return end proc; setXVEC := proc (x::list(
EvalsToFloat)) local i1, D, c; if not nops(x) = N then error 
"SetXVEC expects a list of %1 entries", N end if; if ScalingOn = true then D 
:= Scaling:-getD(); c := Scaling:-getc(); for i1 to N do xvec_swnpsol[i1] := D
[i1]*x[i1]-c[i1] end do else for i1 to N do xvec_swnpsol[i1] := x[i1] end do 
end if; boundsNotViolated(); return end proc; setOption := proc (ReqOpt::(
string = EvalsToFloat)) if not member(lhs(ReqOpt),KnownNpsolOptions) then 
error "requested option %1 not known", lhs(ReqOpt) end if; NpsolOptions[lhs(
ReqOpt)] := rhs(ReqOpt); createOptionsFile(); return end proc; setParameters 
:= proc (NewPars::list(name = EvalsToFloat)) local NewValues, NewParsNames, 
Missing, Obsolete, i1, c, D; NewParsNames := map(lhs,NewPars); Missing, 
Obsolete := Aux:-ListOperations:-getMissingAndObsoleteNames(NewParsNames,
ParNames); if not Missing = {} then error "assignments are missing for %1", 
Missing end if; NewValues := evalf(subs(NewPars,ParNames)); if ScalingParsOn =
true then D := ScalingPars:-getD(); c := ScalingPars:-getc(); for i1 to 
NumParams do Params[i1] := D[i1]*NewValues[i1]-c[i1] end do; return else for 
i1 to NumParams do Params[i1] := NewValues[i1] end do; return end if end proc;
setVariables := proc (vars::list(name = EvalsToFloat)) local GivenVars, 
Missing, Obsolete, InitialPoint; GivenVars := map(lhs,vars); Missing, Obsolete
:= Aux:-ListOperations:-getMissingAndObsoleteNames(GivenVars,RequiredVars); if
not Missing = {} then error "assignments missing for the variables %1", 
Missing end if; InitialPoint := subs(vars,RequiredVars); setXVEC(InitialPoint)
; return end proc; getVarBounds := proc () local LowerBounds, UpperBounds, 
VarBounds, InvD, c, i1; if ScalingOn = true then InvD := Scaling:-getInvD(); c
:= Scaling:-getc(); LowerBounds := array(1 .. N); UpperBounds := array(1 .. N)
; for i1 to N do LowerBounds[i1] := InvD[i1]*(BL_swnpsol[i1]+c[i1]); 
UpperBounds[i1] := InvD[i1]*(BU_swnpsol[i1]+c[i1]) end do; VarBounds := [seq(
RequiredVars[i1] = LowerBounds[i1] .. UpperBounds[i1],i1 = 1 .. N)] else 
VarBounds := [seq(RequiredVars[i1] = BL_swnpsol[i1] .. BU_swnpsol[i1],i1 = 1 
.. N)] end if; return VarBounds end proc; setBoundsForVar := proc (ReqBounds::
{list(name = EvalsToFloat .. EvalsToFloat), name = EvalsToFloat .. 
EvalsToFloat}) local NewBounds, Pos, ValsNewBounds, i1, VarIndex, D, c, 
NewLower, NewUpper, NamesReqBounds, item; if not type(ReqBounds,list) then 
NewBounds := [ReqBounds] else NewBounds := ReqBounds end if; NamesReqBounds :=
map(lhs,ReqBounds); for item in NamesReqBounds while true do if 1 < numboccur(
item,NamesReqBounds) then error 
"symbol %1 occurs more than once in list of bounds", item end if end do; Pos 
:= array(1 .. nops(NewBounds)); for i1 to nops(NewBounds) do try Pos[i1] := 
Aux:-ListOperations:-getPosOfLHSin(lhs(NewBounds[i1]),RequiredVars) catch: 
printf("could not find requested variable %q in list of NLP variables\n",lhs(
NewBounds[i1])); error end try end do; Pos := convert(Pos,list); if ScalingOn
= true then ValsNewBounds := map(rhs,NewBounds); NewLower := [seq(op(1,evalf(
ValsNewBounds[i1])),i1 = 1 .. nops(ValsNewBounds))]; NewUpper := [seq(op(2,
evalf(ValsNewBounds[i1])),i1 = 1 .. nops(ValsNewBounds))]; D := Scaling:-getD(
); c := Scaling:-getc(); for i1 to nops(Pos) do VarIndex := Pos[i1]; 
BL_swnpsol[VarIndex] := D[VarIndex]*NewLower[i1]-c[VarIndex]; BU_swnpsol[
VarIndex] := D[VarIndex]*NewUpper[i1]-c[VarIndex] end do else ValsNewBounds :=
map(rhs,NewBounds); for i1 to nops(Pos) do VarIndex := Pos[i1]; BL_swnpsol[
VarIndex] := op(1,evalf(ValsNewBounds[i1])); BU_swnpsol[VarIndex] := op(2,
evalf(ValsNewBounds[i1])) end do end if; return end proc; 
getIndicesActiveVarBounds := proc () local Eps, ErrsLowerBounds, 
ErrsUpperBounds, i1, OptArgs, ViolatedLower, ViolatedUpper, nLower, nUpper, 
ReqEps; boundsNotViolated(); if nargs < 1 then OptArgs := [] else OptArgs := [
args[1 .. -1]]; if not type(OptArgs,list(name = {EvalsToFloat})) then error 
"valid optional arguments are SetEps= EvalsToFloat " end if; if not `minus`(
map(lhs,convert(OptArgs,set)),{'SetEps'}) = {} then error 
"valid optional arguments are SetEps= EvalsToFloat " end if end if; if 
hasoption(OptArgs,'SetEps',ReqEps,'OptArgs') then Eps := ReqEps else Eps := .1\
e-2 end if; ErrsLowerBounds := array(1 .. N); ErrsUpperBounds := array(1 .. N)
; for i1 to N do ErrsLowerBounds[i1] := abs(2*evalf(xvec_swnpsol[i1]-
BL_swnpsol[i1])/evalf(xvec_swnpsol[i1]+BL_swnpsol[i1])); ErrsUpperBounds[i1] 
:= abs(2*evalf(xvec_swnpsol[i1]-BU_swnpsol[i1])/evalf(xvec_swnpsol[i1]+
BU_swnpsol[i1])) end do; ViolatedLower := array(1 .. N); ViolatedUpper := 
array(1 .. N); nLower := 0; nUpper := 0; for i1 to N do if ErrsLowerBounds[i1]
<= Eps then nLower := nLower+1; ViolatedLower[nLower] := i1 end if; if 
ErrsUpperBounds[i1] <= Eps then nUpper := nUpper+1; ViolatedUpper[nUpper] := 
i1 end if end do; ViolatedLower := convert(ViolatedLower,list); ViolatedUpper
:= convert(ViolatedUpper,list); ViolatedLower := ViolatedLower[1 .. nLower]; 
ViolatedUpper := ViolatedUpper[1 .. nUpper]; return ViolatedLower, 
ViolatedUpper end proc; getIndicesLinearInequalityConstraints := proc () local
i1, nInd, Ind; Ind := array(1 .. NCLIN); nInd := 0; for i1 from N+1 to NCLIN+N
do if not BL_swnpsol[i1] = BU_swnpsol[i1] then nInd := 1+nInd; Ind[nInd] := i1
end if end do; Ind := convert(Ind,list); Ind := Ind[1 .. nInd]; return Ind end
proc; getIndicesNonlinearInequalityConstraints := proc () local i1, nInd, Ind;
Ind := array(1 .. NCNLN); nInd := 0; for i1 from N+NCLIN+1 to N+NCLIN+NCNLN do
if not BL_swnpsol[i1] = BU_swnpsol[i1] then nInd := 1+nInd; Ind[nInd] := i1 
end if end do; Ind := convert(Ind,list); Ind := Ind[1 .. nInd]; return Ind end
proc; getLowerActiveVarBounds := proc () local ViolatedLower, ViolatedUpper, 
VarNames; ViolatedLower, ViolatedUpper := getIndicesActiveVarBounds(args); if
not ViolatedLower = [] then VarNames := [seq(RequiredVars[ViolatedLower[i1]],
i1 = 1 .. nops(ViolatedLower))] else VarNames := [] end if; return VarNames 
end proc; getUpperActiveVarBounds := proc () local ViolatedLower, 
ViolatedUpper, VarNames; ViolatedLower, ViolatedUpper := 
getIndicesActiveVarBounds(args); if not ViolatedUpper = [] then VarNames := [
seq(RequiredVars[ViolatedUpper[i1]],i1 = 1 .. nops(ViolatedUpper))] else 
VarNames := [] end if; return VarNames end proc end module; 
SimpleWrapperTemplate:-createOptionsFile(); return eval(SimpleWrapperTemplate)
end proc; Tests := module () export all, agrawal82, simpleProblem1, 
simpleProblem2; all := proc () if simpleProblem1() = true and simpleProblem1(
adifor) = true and simpleProblem1(adifor,scale) = true and simpleProblem2() =
true and simpleProblem2(adifor) = true and agrawal82() = true and agrawal82(
adifor) = true and agrawal82(adifor,scale) = true then return true else return
false end if end proc; agrawal82 := proc () local anNLP, NPSOLproc, Sol1, 
CalculatedOptimum, ExpectedOptimum, ExpectedObjVal, AcceptedRelErr, i1, RelErr
, Objf; if Aux:-FileOperations:-dirExists("tmpAgrawal82") then error 
"temporary directory ./tmpAgrawal82 needed to run test already exists" else 
mkdir("tmpAgrawal82") end if; anNLP := table(["ExplicitAEs" = [k = 1, K = .12,
a = 5.4, b = 180, c1 = 5, V = 51.3*Pi, mu = k*S*exp(-S/K), sigma = mu/(a+b*S),
c2 = -70.3+270.3*SF], "Parameters" = [F = 3.0, SF = .3], "LinearConstraints" =
[], "Constraints" = [0 = -F/V*X+mu*X, 0 = F/V*(SF-S)-sigma*X], "Variables" = [
X = -infinity .. infinity, S = -infinity .. infinity], "CostFunction" = [c1*X*
F-c2*F]]); anNLP := Aux:-NLP:-parToVarInNLP([F = 0 .. 10, SF = .3 .. 1],anNLP)
; printf("\n"); if 0 < nargs then printf(
"running Agrawal82 with optional argument %s to CreateInstance\n",args[1]); 
NPSOLproc := CreateInstance(anNLP,"./tmpAgrawal82",args[1 .. -1]) else printf(
"running Agrawal82 without option adifor\n"); NPSOLproc := CreateInstance(
anNLP,"./tmpAgrawal82") end if; Sol1 := [F = 3, SF = .3, X = 2.6201, S = .\
22443e-1]; NPSOLproc:-setVariables(Sol1); NPSOLproc:-runOpt(); Objf := 
NPSOLproc:-getObjf(); CalculatedOptimum := NPSOLproc:-getVariables(); dlclose(
"cwrap_npsol1.so"); system("rm -r tmpAgrawal82"); ExpectedOptimum := [X = 16.4\
4085444, S = .6800249428e-1, F = 6.218460030, SF = .9999999997]; 
ExpectedObjVal := -732.508024836933600; AcceptedRelErr := .1e-2; for i1 to 
nops(ExpectedOptimum) do if not abs(2*evalf(rhs(CalculatedOptimum[i1])-rhs(
ExpectedOptimum[i1]))/evalf(rhs(CalculatedOptimum[i1])+rhs(ExpectedOptimum[i1]
))) < AcceptedRelErr then return false end if end do; RelErr := abs(2*evalf(
ExpectedObjVal-Objf)/evalf(ExpectedObjVal+Objf)); if not RelErr < 
AcceptedRelErr then return false end if; return true end proc; simpleProblem1
:= proc () local g, xRange, yRange, anNLP, StandardNLP, NPSOLproc, Objf, 
OptPoint, ExpectedObjf, ExpectedOptPoint, Inform, i1, RelErr, rmv, argslist; 
argslist := [args]; if member(noremove,argslist) then rmv := false; Aux:-
ListOperations:-removeItemFromList(noremove,argslist) else rmv := true end if;
ExpectedObjf := 1.00000; ExpectedOptPoint := [1.0, 1.0]; if Aux:-
FileOperations:-dirExists("tmpSimpleProblem1") then error 
"temporary directory ./tmpSimpleProblem1 needed to run test already exists" 
else mkdir("tmpSimpleProblem1") end if; g := (x-2)^2+(y-1)^2; xRange := -
infinity .. infinity; yRange := -infinity .. infinity; anNLP := Aux:-NLP:-
newNLP(); anNLP["CostFunction"] := [eval(g)]; anNLP["Variables"] := [x = 
xRange, y = yRange]; anNLP["Constraints"] := [0 < y-x^2]; anNLP[
"LinearConstraints"] := [0 < 2-x-y]; anNLP["ExplicitAEs"] := []; printf("\n");
printf("running SimpleProblem1\n"); if 0 < nargs then NPSOLproc := NPSOL:-
CreateInstance(anNLP,"./tmpSimpleProblem1",op(argslist[1 .. -1])) else 
NPSOLproc := NPSOL:-CreateInstance(anNLP,"./tmpSimpleProblem1") end if; 
NPSOLproc:-setXVEC([0, 1]); NPSOLproc:-runOpt(); Objf := NPSOLproc:-getObjf();
OptPoint := NPSOLproc:-getXVEC(); Inform := NPSOLproc:-getInform(); for i1 to
nops(anNLP["Variables"]) do RelErr := abs(2*evalf(ExpectedOptPoint[i1]-
OptPoint[i1])/evalf(ExpectedOptPoint[i1]+OptPoint[i1])); if not RelErr < .1e-4
then return false end if end do; RelErr := abs(2*evalf(ExpectedObjf-Objf)/
evalf(ExpectedObjf+Objf)); if not RelErr < .1e-2 then return false end if; 
unassign('NPSOLproc'); NPSOLproc := RestoreInstance("./tmpSimpleProblem1"); 
NPSOLproc:-setXVEC([0, 1]); NPSOLproc:-runOpt(); Objf := NPSOLproc:-getObjf();
OptPoint := NPSOLproc:-getXVEC(); Inform := NPSOLproc:-getInform(); if rmv 
then dlclose("cwrap_npsol1.so"); system("rm -r tmpSimpleProblem1") end if; for
i1 to nops(anNLP["Variables"]) do RelErr := abs(2*evalf(ExpectedOptPoint[i1]-
OptPoint[i1])/evalf(ExpectedOptPoint[i1]+OptPoint[i1])); if not RelErr < .1e-4
then return false end if end do; RelErr := abs(2*evalf(ExpectedObjf-Objf)/
evalf(ExpectedObjf+Objf)); if not RelErr < .1e-2 then return false end if; 
return true end proc; simpleProblem2 := proc () local g, xRange, yRange, anNLP
, StandardNLP, NPSOLproc, Objf, OptPoint, ExpectedObjf, ExpectedOptPoint, 
Inform, RelErr, i1, AbsErr; ExpectedObjf := 0.; ExpectedOptPoint := [1.0, 1.0]
; if Aux:-FileOperations:-dirExists("tmpSimpleProblem2") then error 
"temporary directory ./tmpSimpleProblem2 needed to run test already exists" 
else mkdir("tmpSimpleProblem2") end if; CreateStandardOptionsFile(
"tmpSimpleProblem2"); g := 100*(y-x^2)^2+(1-x)^2; xRange := -infinity .. 
infinity; yRange := -infinity .. infinity; anNLP := Aux:-NLP:-newNLP(); anNLP[
"CostFunction"] := [eval(g)]; anNLP["Variables"] := [x = xRange, y = yRange];
anNLP["ExplicitAEs"] := []; printf("\n"); printf("running SimpleProblem2\n");
if 0 < nargs then NPSOLproc := CreateInstance(anNLP,"./tmpSimpleProblem2",args
[1 .. -1]) else printf("running SimpleProblem2 without option adifor\n"); 
NPSOLproc := CreateInstance(anNLP,"./tmpSimpleProblem2") end if; NPSOLproc:-
setXVEC([-1.2, 1.0]); NPSOLproc:-runOpt(); Objf := NPSOLproc:-getObjf(); 
OptPoint := NPSOLproc:-getXVEC(); Inform := NPSOLproc:-getInform(); for i1 to
nops(anNLP["Variables"]) do RelErr := abs(2*evalf(ExpectedOptPoint[i1]-
OptPoint[i1])/evalf(ExpectedOptPoint[i1]+OptPoint[i1])); if not RelErr < .1e-2
then return false end if end do; AbsErr := abs(ExpectedObjf-Objf); if not 
AbsErr < .1e-2 then return false end if; NPSOLproc:-setXVEC([-.3, .1]); 
NPSOLproc:-runOpt(); Objf := NPSOLproc:-getObjf(); OptPoint := NPSOLproc:-
getXVEC(); Inform := NPSOLproc:-getInform(); dlclose("cwrap_npsol1.so"); 
system("rm -r tmpSimpleProblem2"); for i1 to nops(anNLP["Variables"]) do 
RelErr := abs(2*evalf(ExpectedOptPoint[i1]-OptPoint[i1])/evalf(
ExpectedOptPoint[i1]+OptPoint[i1])); if not RelErr < .1e-2 then return false 
end if end do; RelErr := abs(ExpectedObjf-Objf); if not AbsErr < .1e-2 then 
return false end if; return true end proc; end module; end module;
