BVPSol := module () export createExternalFunctions, CreateInstance, 
createSharedObject, runDefineExternal, SubsIntoSimpleWrapper, Tests; 
createExternalFunctions := proc (anBVPSys::DAESys) local n, npars, rhsODEs, 
FCN, listOfSubs, i, rhsOfOdesWithSubs, st1, st2, st3, st4, stFCN, fdFCN, 
fdFCN2, filePathFCN, filePathFCNFortran, stDFDY, j, jac, fdDFDY, filePathDFDY,
filePathDFDYFortran, jacWithSubs, DFDY; anBVPSys := Aux:-SystemClasses:-
subsExplicitAEsIntoDAESys(anBVPSys); n := nops(anBVPSys["ODEs"]); npars := 
nops(anBVPSys["Parameters"]); rhsODEs := map(rhs,anBVPSys["ODEs"]); listOfSubs
:= []; for i to n do listOfSubs := [op(listOfSubs), anBVPSys["DynVars"][i] = Z
(i)] end do; for i to npars do listOfSubs := [op(listOfSubs), lhs(anBVPSys[
"Parameters"][i]) = TRPAR(i)] end do; rhsOfOdesWithSubs := subs(listOfSubs,
rhsODEs); st1 := "proc(T,Z,DZ,TRPAR) \n"; st2 := StringTools[Join]([
" Z:=Vector(", convert(n,string), ");\n DZ:=Vector(", convert(n,string), 
");\n TRPAR:=Vector(", convert(npars,string), ");\n"],""); st3 := ""; for i to
n do st3 := StringTools[Join]([st3, " DZ(", convert(i,string), "):=", convert(
rhsOfOdesWithSubs[i],string), ";\n"],"") end do; st4 := " return;\nend proc;";
stFCN := StringTools[Join]([st1, st2, st3, st4],""); fdFCN := fopen(
CFileForFCN,WRITE); fprintf(fdFCN,stFCN); fclose(fdFCN); filePathFCN := cat(
currentdir(),"/CFileForFCN"); if not Aux:-FileOperations:-fileExists(
filePathFCN) then error "file with C code for FCN procedure was not created" 
else filePathFCNFortran := cat(currentdir(),"/fortranrhs.f"); if Aux:-
FileOperations:-fileExists(filePathFCNFortran) then system(cat(`rm -f  `,
"fortranrhs.f")) end if; read filePathFCN; FCN := %; CodeGeneration[Fortran](
FCN,output = "fortranrhs.f"); system(cat(`rm -f  `,"CFileForFCN")) end if; jac
:= VectorCalculus[Jacobian](rhsODEs,anBVPSys["DynVars"]); jacWithSubs := subs(
listOfSubs,jac); st1 := "proc(T,Z,DF,TRPAR) \n"; st2 := StringTools[Join]([
" Z:=Vector(", convert(n,string), ");\n DF:=Matrix(", convert(n,string), ",",
convert(n,string), ");\n TRPAR:=Vector(", convert(npars,string), ");\n"],"");
st3 := ""; for i to n do for j to n do st3 := StringTools[Join]([st3, " DF(",
convert(i,string), ",", convert(j,string), "):=", convert(jacWithSubs[i,j],
string), ";\n"],"") end do end do; stDFDY := StringTools[Join]([st1, st2, st3,
st4],""); fdDFDY := fopen(CFileForDFDY,WRITE); fprintf(fdDFDY,stDFDY); fclose(
fdDFDY); filePathDFDY := cat(currentdir(),"/CFileForDFDY"); if not Aux:-
FileOperations:-fileExists(filePathDFDY) then error 
"file with C code for DFDY procedure was not created" else filePathDFDYFortran
:= cat(currentdir(),"/fortranjac.f"); if Aux:-FileOperations:-fileExists(
filePathDFDYFortran) then system(cat(`rm -f  `,"fortranjac.f")) end if; read 
filePathDFDY; DFDY := %; CodeGeneration[Fortran](DFDY,output = "fortranjac.f")
; system(cat(`rm -f  `,"CFileForDFDY")) end if; RETURN() end proc; 
CreateInstance := proc (anBVPSys::DAESys, WorkingDir::string, nodes::integer,
objNames::list(string)) local bvpsol, BVPSolInstance; if Aux:-FileOperations:-
fileExists(cat(WorkingDir,"/cwrap_bvpsol1.so")) then error "in working dir, %1\
, a maple wrapper file named cwrap_bvpsol1.so already exists", WorkingDir end
if; anBVPSys := Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(anBVPSys); try 
bvpsol := createSharedObject(WorkingDir,anBVPSys,objNames) catch: printf(
"CreateSharedObject failed, \n"); printf("  error thrown was %q\n",
lastexception); printf("  attempting to continue CreateInstance\n") end try; 
printf("creating instance of interface\n"); try BVPSolInstance := 
SubsIntoSimpleWrapper(anBVPSys,bvpsol,WorkingDir,nodes) catch: printf(
"SubsIntoSimpleWrapper failed, \n"); printf("  error thrown was %q\n",
lastexception); printf("  attempting to continue CreateInstance\n") end try; 
RETURN(eval(BVPSolInstance)) end proc; createSharedObject := proc (bvpsolDir::
string, anBVPSys::DAESys, objNames::list(string)) local oldDir, LinkerCommand,
bvpsol, SharedObjectName, InstanceCounter, i, objNamesList; InstanceCounter :=
1; if not Aux:-FileOperations:-dirExists(bvpsolDir) then mkdir(bvpsolDir) end
if; oldDir := currentdir(); currentdir(bvpsolDir); createExternalFunctions(
anBVPSys); system("gfortran -c -shared -O -fPIC fortranrhs.f -o fortranrhs.o")
; system("gfortran -c -shared -O -fPIC fortranjac.f -o fortranjac.o"); 
SharedObjectName := cat("cwrap_bvpsol",InstanceCounter,".so"); objNamesList :=
""; for i to nops(objNames) do objNamesList := cat(objNamesList,objNames[i],
" ") end do; LinkerCommand := cat(
`gcc -shared -Xlinker -Bsymbolic fortranrhs.o fortranjac.o `,_ModulesDirectory
,`/BVPSol/ext_routines/f2c/cwrap/bvpopen.o `,_ModulesDirectory,
`/BVPSol/ext_routines/f2c/cwrap/bvpclose.o `,_ModulesDirectory,
`/BVPSol/ext_routines/f2c/cwrap/cwrap_bvpsol.o `,objNamesList,
_ModulesDirectory,`/BVPSol/ext_routines/f2c/shared_obj/period.so `,
` -lf2c -O2 -lm -lgfortran -lTIDES -o `,SharedObjectName); system(
LinkerCommand); bvpsol := runDefineExternal(SharedObjectName); currentdir(
oldDir); RETURN(bvpsol) end proc; runDefineExternal := proc (SharedObjectName
::string) local bvpsol; bvpsol := define_external('cwrapper_bvpsol',N::integer
[4],M::integer[4],IPRINT::integer[4],X0::ARRAY(1 .. N,float[8]),P::float[8],
NPAR::integer[4],TRPAR::ARRAY(1 .. NPAR,float[8]),Y::ARRAY(1 .. M,1 .. N,float
[8]),FM::ARRAY(1 .. 2,1 .. N,float[8]),POUT::REF(float[8]),IFAIL::REF(integer[
4]),ERRY::REF(float[8]),FAL::ARRAY(1 .. N,1 .. NPAR,float[8],order = 
Fortran_order),FXAL::ARRAY(1 .. N,1 .. N,1 .. NPAR,float[8],order = 
Fortran_order),FXX::ARRAY(1 .. N,1 .. N,1 .. N,float[8],order = Fortran_order)
,FX::ARRAY(1 .. N,1 .. N,float[8],order = Fortran_order),HDIF::float[8],LIB =
SharedObjectName); return bvpsol end proc; SubsIntoSimpleWrapper := proc (
anBVPSys::DAESys, bvpsolCwrapper::procedure, WorkingDir::string, nodes::
integer) local N, M, IPRINT, NPAR, X0, TRPAR, i, P, Y, FM, POUT, IFAIL, ERRY,
FAL, FXAL, FXX, listForFXAL, listForFXX, SimpleWrapperTemplate, FX, HDIF, 
FAlAl, FXXAl, FXAlAl, listForFXXPdim1, listForFXXPdim2, listForFXPPdim1, 
listForFXPPdim2; N := nops(anBVPSys["ODEs"]); M := nodes; NPAR := nops(
anBVPSys["Parameters"]); IPRINT := 0; X0 := Vector(N,datatype = float[8],order
= C_order); TRPAR := Vector(NPAR,datatype = float[8],order = C_order); P := 0;
Y := Matrix(M,N,datatype = float[8],order = C_order); FM := Matrix(2,N,
datatype = float[8],order = C_order); FAL := Matrix(N,NPAR,datatype = float[8]
,order = Fortran_order); listForFXAL := []; for i to N do listForFXAL := [op(
listForFXAL), convert(Matrix(N,NPAR),listlist)] end do; FXAL := Array(
listForFXAL,datatype = float[8],order = Fortran_order); listForFXX := []; for
i to N do listForFXX := [op(listForFXX), convert(Matrix(N,N),listlist)] end do
; FXX := Array(listForFXX,datatype = float[8],order = Fortran_order); FX := 
Matrix(N,N,datatype = float[8],order = Fortran_order); HDIF := .1e-2; 
SimpleWrapperTemplate := module () local theBVPSys; export getEigVals, getERRY
, getFalpha, getFxalpha, getFx, getFxx, getIFAIL, getInitPeriod, getInitPoint,
getIPRINT, getM, getN, getNPAR, getParams, getPeriod, getSys, getSol, 
getStepSize, getXT, runBVPSol, setInitPeriod, setInitPoint, setIPRINT, 
setParams, setStepSize; theBVPSys := copy(anBVPSys); getEigVals := proc () 
return eval(FM) end proc; getERRY := proc () return eval(ERRY) end proc; 
getFalpha := proc () return eval(FAL) end proc; getFxalpha := proc () return 
eval(FXAL) end proc; getFx := proc () return eval(FX) end proc; getFxx := proc
() return eval(FXX) end proc; getIFAIL := proc () return eval(IFAIL) end proc;
getInitPeriod := proc () return eval(P) end proc; getInitPoint := proc () 
return eval(X0) end proc; getIPRINT := proc () return eval(IPRINT) end proc; 
getM := proc () return eval(M) end proc; getN := proc () return eval(N) end 
proc; getNPAR := proc () return eval(NPAR) end proc; getParams := proc () 
return eval(TRPAR) end proc; getPeriod := proc () return eval(POUT) end proc;
getSys := proc () return eval(theBVPSys) end proc; getSol := proc () return 
eval(Y) end proc; getStepSize := proc () return eval(HDIF) end proc; getXT :=
proc () local i, XT; if 0 <= POUT and 0 < M then XT := Vector(M); for i to M 
do XT[i] := POUT*(i-1)/M end do end if; return eval(XT) end proc; runBVPSol :=
proc () local NextNumber, NextFileName, OldDir; OldDir := currentdir(); 
currentdir(WorkingDir); bvpsolCwrapper(N,M,IPRINT,X0,P,NPAR,TRPAR,Y,FM,'POUT',
'IFAIL','ERRY',FAL,FXAL,FXX,FX,HDIF); currentdir(OldDir); return end proc; 
setInitPeriod := proc (periodGuess::EvalsToFloat) P := periodGuess; return end
proc; setInitPoint := proc (NewPars::list(name = EvalsToFloat)) local 
NewValues, NewParsNames, Missing, Obsolete, i1, ParNames; ParNames := 
theBVPSys["DynVars"]; NewParsNames := map(lhs,NewPars); Missing, Obsolete := 
Aux:-ListOperations:-getMissingAndObsoleteNames(NewParsNames,ParNames); if not
Missing = {} then error "assignments are missing for %1", Missing end if; 
NewValues := subs(NewPars,ParNames); for i1 to nops(ParNames) do X0[i1] := 
NewValues[i1] end do; return end proc; setIPRINT := proc (valIPRINT::integer)
if not (valIPRINT = -1 or valIPRINT = 0 or valIPRINT = 1) then error 
"IRINT can be set only to -1, 0, or 1" end if; IPRINT := valIPRINT; return end
proc; setParams := proc (NewPars::list(name = EvalsToFloat)) local NewValues,
NewParsNames, Missing, Obsolete, i1, ParNames; ParNames := map(lhs,theBVPSys[
"Parameters"]); NewParsNames := map(lhs,NewPars); Missing, Obsolete := Aux:-
ListOperations:-getMissingAndObsoleteNames(NewParsNames,ParNames); if not 
Missing = {} then error "assignments are missing for %1", Missing end if; 
NewValues := subs(NewPars,ParNames); for i1 to nops(ParNames) do TRPAR[i1] :=
NewValues[i1] end do; return end proc; setStepSize := proc (stepH::
EvalsToFloat) HDIF := stepH; return end proc end module; return eval(
SimpleWrapperTemplate) end proc; Tests := module () export problemOlsen1983; 
problemOlsen1983 := proc () local expectedPeriod, aSys, BVPSOLproc, eigVals, 
period, nodes, sol, succeeded, listOfTIDESObj; expectedPeriod := 14.5391926329\
389083; if Aux:-FileOperations:-dirExists("tmpOlsen") then error 
"temporary directory ./tmpOlsen needed to run test already exists" else mkdir(
"tmp") end if; aSys := Models:-Olsen1983:-Sys:-getSys(); listOfTIDESObj := [
cat(_ModulesDirectory,"/BVPSol/demo/ext_routines/fromTIDES/olsen/olsen.o"), 
cat(_ModulesDirectory,"/BVPSol/demo/ext_routines/fromTIDES/olsen/calltides.o")
]; BVPSOLproc := BVPSol:-CreateInstance(aSys,"./tmpOlsen",30,listOfTIDESObj);
BVPSOLproc:-setParams([k1 = .19, k3 = .66651056e-2]); BVPSOLproc:-setInitPoint
([A = 6.0, B = 195.0, X = .3e-1, Y = .7e-1]); BVPSOLproc:-setInitPeriod(14.0);
BVPSOLproc:-runBVPSol(); period := BVPSOLproc:-getPeriod(); eigVals := 
BVPSOLproc:-getEigVals(); nodes := BVPSOLproc:-getXT(); sol := BVPSOLproc:-
getSol(); succeeded := BVPSOLproc:-getIFAIL(); if not ((period-expectedPeriod)
^2 < 1/10000000 and 0 < succeeded) then printf(
"The test problem is not succeeded"); return false end if; printf(
"The test problem is succeeded"); return true end proc; end module; end module
;
