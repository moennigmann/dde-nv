Aux := module () export BoxOperations, CreateCurve, Derivs, FileOperations, 
equalWithinRelErr, equalWithinFudgedRelErr, getAbsErr, getFudgedRelErr, 
getRelErr, init, IntervalArithmetics, LinearAlgebra, LinearEqns, 
ListOperations, NLP, Other, Parsers, Programming, SystemClasses, TensProd, 
TransferNlpToGams; global `type/EvalsToFloat`, `type/NLP`, `type/BCNLP`, 
`type/term`, `type/DTASys`, `type/DAESys`; option package; BoxOperations := 
module () export getIndexToLongestAxis; getIndexToLongestAxis := proc (box::
list([numeric, numeric])) local i, width, spaceDimension, indexToMax, maxWidth
; spaceDimension := nops(box); maxWidth := box[1][2]-box[1][1]; indexToMax :=
1; for i from 2 to spaceDimension do width := box[i][2]-box[i][1]; if maxWidth
< width then maxWidth := width; indexToMax := i end if end do; return 
indexToMax end proc; end module; CreateCurve := proc (PointLikeData::table, 
OtherData::list(name = anything)) local item, MandatoryEntries, 
ExistingEntries, NumDataPointsInConstr, VarNamesInConstr, SubsListParsInConstr
, TemplateModule; MandatoryEntries := [Parameters, NumPoints]; ExistingEntries
:= map(lhs,OtherData); for item in MandatoryEntries while true do if not 
member(item,ExistingEntries) then error "second argument must have entry %1",
item end if end do; SubsListParsInConstr := Aux:-ListOperations:-getRHSofIn(
Parameters,OtherData); if not type(SubsListParsInConstr,{list(name = 
EvalsToFloat), []}) then error "entry Parameters in second argument is expecte\
d to be of type {[], list(name= EvalsToFloat)}" end if; NumDataPointsInConstr
:= Aux:-ListOperations:-getRHSofIn('NumPoints',OtherData); VarNamesInConstr :=
map(op,[indices(PointLikeData)]); for item in VarNamesInConstr while true do 
if not type(PointLikeData[item],array) then error "all entries in 1st argument\
 except Parameters and NumPoints must be arrays, check failed for %1", item 
end if end do; TemplateModule := module () local Indices, NumDataPoints, 
VarNames, data, ParNames, getSinglePlot, SubsListPars, ExtraData; export 
createPS, extendCurveByEAEs, getData, getNonTableData, getNumPoints, 
getParameters, getParNames, getSinglePlotForFunction, 
getSinglePlotForFunction3d, getVariables, getVarNames, getPoint, getPlot, 
joinCurve, splitCurve; VarNames := VarNamesInConstr; NumDataPoints := 
NumDataPointsInConstr; data := copy(PointLikeData); ExtraData := OtherData; 
SubsListPars := SubsListParsInConstr; ParNames := map(lhs,Aux:-ListOperations
:-getRHSofIn(Parameters,ExtraData)); getData := proc () return eval(
PointLikeData) end proc; getParameters := proc () local Pars; Pars := Aux:-
ListOperations:-getRHSofIn(Parameters,ExtraData); return Pars end proc; 
getPlot := proc (xName::name, yName::{name, term, list(name)}) local Result, 
i1, PlotOpts, LegendString, xNameString; xNameString := convert(xName,string);
PlotOpts := [style = line]; if type(yName,list(name)) then LegendString := map
(convert,yName,string); Result := [seq(getSinglePlot(xName,yName[i1],op(
PlotOpts),linestyle = i1,labels = [xNameString, ""],legend = LegendString[i1],
args[3 .. -1]),i1 = 1 .. nops(yName))] elif type(yName,name) then LegendString
:= convert(yName,string); Result := getSinglePlot(xName,yName,op(PlotOpts),
linestyle = 1,args[3 .. -1]) elif type(yName,term) then Result := 
getSinglePlotForFunction(xName,yName,op(PlotOpts),linestyle = 1,args[3 .. -1])
end if; return Result end proc; getSinglePlot := proc (xName::name, yName::
name) local points, ListOfPoints, ListOfPlots, PlotObj; if not member(xName,{
`#`, op(VarNames)}) then error "%1 must be a variable or #", xName end if; if
not member(yName,VarNames) then error "%1 is not a variable", yName end if; if
not xName = `#` then points := [seq([PointLikeData[xName][i1], PointLikeData[
yName][i1]],i1 = 0 .. NumDataPoints)] else points := [seq([i1, PointLikeData[
yName][i1]],i1 = 0 .. NumDataPoints)] end if; PlotObj := plots[pointplot](
points,args[3 .. -1]); return PlotObj end proc; getSinglePlotForFunction := 
proc (xName::name, yName::term) local points, ListOfPoints, ListOfPlots, 
PlotObj, ListOfIndets, item, AllNames, i1, i2, NewPoint, NumIndets, 
VarsInIndets, ParsInIndets, NumVarsInIndets, SubsListPars, SubsList, xValue; 
if not member(xName,{`#`, op(VarNames)}) then error 
"%1 must be a variable or #", xName end if; ListOfIndets := Aux:-
ListOperations:-getListOfIndetsIn(yName); NumIndets := nops(ListOfIndets); 
AllNames := [op(VarNames), op(ParNames)]; for item in ListOfIndets while true
do if not member(item,AllNames) then error 
"%1 is not a variable, nor a parameter ", item end if end do; VarsInIndets :=
convert(`intersect`(convert(VarNames,set),convert(ListOfIndets,set)),list); 
ParsInIndets := convert(`intersect`(convert(ParNames,set),convert(ListOfIndets
,set)),list); NumVarsInIndets := nops(VarsInIndets); SubsListPars := subs(
GetParameters(),ParsInIndets); SubsListPars := [seq(ParsInIndets[i1] = 
SubsListPars[i1],i1 = 1 .. nops(ParsInIndets))]; points := array(1 .. 
NumDataPoints); for i1 to NumDataPoints do SubsList := [seq(VarsInIndets[i2] =
PointLikeData[VarsInIndets[i2]][i1],i2 = 1 .. NumVarsInIndets)]; if not xName
= `#` then xValue := PointLikeData[xName][i1] else xValue := i1 end if; 
NewPoint := [xValue, evalf(subs(SubsList,SubsListPars,yName))]; points[i1] :=
NewPoint end do; points := convert(points,list); PlotObj := plots[pointplot](
points,args[3 .. -1]); return PlotObj end proc; getSinglePlotForFunction3d :=
proc (xName::name, yName::{name, numeric}, func::term) local points, 
ListOfPoints, ListOfPlots, PlotObj, ListOfIndets, item, AllNames, i1, i2, 
NewPoint, NumIndets, VarsInIndets, ParsInIndets, NumVarsInIndets, SubsListPars
, SubsList, xValue, yValue; if not member(xName,{`#`, op(VarNames)}) then 
error "%1 must be a variable or #", xName end if; if not (member(yName,{`#`, 
op(VarNames)}) or type(yName,numeric)) then error "%1 must be a variable or #"
, yName end if; ListOfIndets := Aux:-ListOperations:-getListOfIndetsIn(func);
NumIndets := nops(ListOfIndets); AllNames := [op(VarNames), op(ParNames)]; for
item in ListOfIndets while true do if not member(item,AllNames) then error 
"%1 is not a variable, nor a parameter ", item end if end do; VarsInIndets :=
convert(`intersect`(convert(VarNames,set),convert(ListOfIndets,set)),list); 
ParsInIndets := convert(`intersect`(convert(ParNames,set),convert(ListOfIndets
,set)),list); NumVarsInIndets := nops(VarsInIndets); SubsListPars := subs(
GetParameters(),ParsInIndets); SubsListPars := [seq(ParsInIndets[i1] = 
SubsListPars[i1],i1 = 1 .. nops(ParsInIndets))]; points := array(1 .. 
NumDataPoints); for i1 to NumDataPoints do SubsList := [seq(VarsInIndets[i2] =
PointLikeData[VarsInIndets[i2]][i1],i2 = 1 .. NumVarsInIndets)]; if not xName
= `#` then xValue := PointLikeData[xName][i1] else xValue := i1 end if; if not
yName = `#` then yValue := PointLikeData[yName][i1] else yValue := i1 end if;
NewPoint := [xValue, yValue, evalf(subs(SubsList,SubsListPars,func))]; points[
i1] := NewPoint end do; points := convert(points,list); PlotObj := plots[
pointplot3d](points,args[4 .. -1]); return PlotObj end proc; createPS := proc
(xName::name, yName::{name, list(name)}, FileName::string) local p; if Aux:-
FileOperations:-fileExists(FileName) then error "file %1 already exists", 
FileName end if; p := getPlot(args[1 .. 2],args[4 .. nargs]); plotsetup(ps,
plotoutput = FileName); print(p); plotsetup(default); return end proc; 
getNonTableData := proc () return ExtraData end proc; getNumPoints := proc ()
return NumDataPoints end proc; getPoint := proc (ReqNum::posint) local i1, 
point; if not ReqNum <= NumDataPoints then error 
"only %1 points are available", NumDataPoints end if; point := [op(
getParameters()), seq(VarNames[i1] = data[VarNames[i1]][ReqNum],i1 = 1 .. nops
(VarNames))]; return point end proc; getVarNames := proc () return VarNames 
end proc; getParNames := proc () return ParNames end proc; joinCurve := proc (
curve2::`module`) local SetNewVars, SetOldVars, SetNewPars, SetOldPars, 
ObsoleteInNewCurve, ObsoleteInOldCurve, NewNumPoints, NewDataList, DataCurve2,
NewDataTable, item, i1, NewCurve; SetNewVars := curve2:-getVarNames(); 
SetNewVars := convert(SetNewVars,set); SetOldVars := convert(VarNames,set); 
ObsoleteInNewCurve := `minus`(SetNewVars,SetOldVars); if not 
ObsoleteInNewCurve = {} then error "first argument has obsolete variables %1",
ObsoleteInNewCurve end if; ObsoleteInOldCurve := `minus`(SetOldVars,SetNewVars
); if not ObsoleteInOldCurve = {} then error 
"first argument is missing variables %1", ObsoleteInOldCurve end if; 
SetNewPars := curve2:-getParNames(); SetNewPars := convert(SetNewPars,set); 
SetOldPars := convert(ParNames,set); ObsoleteInNewCurve := `minus`(SetNewPars,
SetOldPars); if not ObsoleteInNewCurve = {} then error 
"first argument has obsolete parameters %1", ObsoleteInNewCurve end if; 
ObsoleteInOldCurve := `minus`(SetOldPars,SetNewPars); if not 
ObsoleteInOldCurve = {} then error "first argument is missing parameters %1",
ObsoleteInOldCurve end if; if not Aux:-ListOperations:-noConflictsInSubsLists(
getParameters(),curve2:-getParameters()) then error cat(
"curves cannot be joined since there are ","conflicts w.r.t. Parameters") end
if; NewNumPoints := NumDataPoints+curve2:-getNumPoints(); NewDataList := 
ExtraData; NewDataList := Aux:-ListOperations:-setRHSofInTo(NumPoints,
NewDataList,NewNumPoints); DataCurve2 := curve2:-getData(); NewDataTable := 
table(); for item in VarNames while true do NewDataTable[item] := array(1 .. 
NewNumPoints); for i1 to NumDataPoints do NewDataTable[item][i1] := data[item]
[i1] end do; for i1 to curve2:-getNumPoints() do NewDataTable[item][
NumDataPoints+i1] := DataCurve2[item][i1] end do end do; NewCurve := 
CreateCurve(NewDataTable,NewDataList); return eval(NewCurve) end proc; 
splitCurve := proc (DecisionMaker::procedure) local i1, IndicesFutureCurve1, 
IndicesFutureCurve2, NewCurve1, NewCurve2, NumPointsCurve1, NumPointsCurve2, 
DataCurve1, DataCurve2, item, ListCurve1, ListCurve2; IndicesFutureCurve1 := [
]; IndicesFutureCurve2 := []; for i1 to NumDataPoints do if DecisionMaker(i1,
getPoint(i1)) = true then IndicesFutureCurve1 := [op(IndicesFutureCurve1), i1]
else IndicesFutureCurve2 := [op(IndicesFutureCurve2), i1] end if end do; if 
IndicesFutureCurve1 = [] then NewCurve1 := Aux:-CreateCurve(data,ExtraData); 
return [eval(NewCurve1), 'NewCurve2'] elif IndicesFutureCurve2 = [] then 
NewCurve2 := Aux:-CreateCurve(data,ExtraData); return ['NewCurve1', eval(
NewCurve2)] end if; NumPointsCurve1 := nops(IndicesFutureCurve1); 
NumPointsCurve2 := nops(IndicesFutureCurve2); DataCurve1 := table(); for item
in VarNames while true do DataCurve1[item] := array(1 .. NumPointsCurve1); for
i1 to NumPointsCurve1 do DataCurve1[item][i1] := data[item][
IndicesFutureCurve1[i1]] end do end do; ListCurve1 := ExtraData; ListCurve1 :=
Aux:-ListOperations:-setRHSofInTo(NumPoints,ListCurve1,NumPointsCurve1); 
DataCurve2 := table(); for item in VarNames while true do DataCurve2[item] :=
array(1 .. NumPointsCurve2); for i1 to NumPointsCurve2 do DataCurve2[item][i1]
:= data[item][IndicesFutureCurve2[i1]] end do end do; ListCurve2 := ExtraData;
ListCurve2 := Aux:-ListOperations:-setRHSofInTo(NumPoints,ListCurve2,
NumPointsCurve2); NewCurve1 := Aux:-CreateCurve(DataCurve1,ListCurve1); 
NewCurve2 := Aux:-CreateCurve(DataCurve2,ListCurve2); return [NewCurve1, 
NewCurve2] end proc; getVariables := proc (ReqNum::posint) local i1, point; if
not ReqNum <= NumDataPoints then error "only %1 points are available", 
NumDataPoints end if; point := [seq(VarNames[i1] = data[VarNames[i1]][ReqNum],
i1 = 1 .. nops(VarNames))]; return point end proc; extendCurveByEAEs := proc (
Sys::{DAESys, ExtAESys, NLP}, WorkingDir::string) local ListOfErrs, LhsEqns, 
NumEAEs, ErrPars, ErrVars, EAEsNamesInEAEs, NumEAEsNamesInEAEs, NumVarNames, 
item, ValsPars, i1, ValsVars, AllVals, ValsEAEs, NewData, NewCurve, Eqns; Eqns
:= Sys[ExplicitAEs]; LhsEqns := convert(map(lhs,Eqns),set); ErrPars := 
`intersect`(LhsEqns,convert(ParNames,set)); if not ErrPars = {} then error 
"equations given contain lhs %1 which is parameter of curve", ErrPars end if;
ErrVars := `intersect`(LhsEqns,convert(ParNames,set)); if not ErrVars = {} 
then error "equations given contain lhs %1 which is variable of curve", 
ErrVars end if; if not type(NumSys,`module`) then error 
"module NumSys is not available " end if; NumEAEs := NumSys:-EAE:-
CreateInstance(Sys,WorkingDir,'new','BuildInterface'); EAEsNamesInEAEs := 
NumEAEs:-GetEAEsNames(); NumEAEsNamesInEAEs := nops(EAEsNamesInEAEs); 
NumVarNames := nops(VarNames); NewData := copy(data); for item in 
EAEsNamesInEAEs while true do NewData[item] := array(1 .. NumDataPoints) end 
do; ValsPars := GetParameters(); for i1 to NumDataPoints do ValsVars := 
GetVariables(i1); AllVals := [op(ValsVars), op(ValsPars)]; NumEAEs:-
SetVariables(AllVals); NumEAEs:-SetParameters(AllVals); ValsEAEs := NumEAEs:-
GetExplicitAEs(); for item in EAEsNamesInEAEs while true do NewData[item][i1]
:= subs(ValsEAEs,item) end do end do; NewCurve := Aux:-CreateCurve(NewData,
GetNonTableData()); return NewCurve end proc end module; return eval(
TemplateModule) end proc; Derivs := module () export f_p, f_x, f_xx, f_xp; f_p
:= proc (functions::{Vector, list}, p::{Vector, list}) local f, Result; if 
type(functions,list(equation)) then f := map(rhs,functions) else f := 
functions end if; Result := VectorCalculus[Jacobian](f,p); return eval(Result)
end proc; f_x := proc (functions::{Vector, list}, x::{Vector, list}) local nx,
nf, Result; nx := nops(x); nf := nops(functions); if not nf = nx then error 
"dimensions of f and x are not the same" end if; Result := Derivs:-f_p(
functions,x); return eval(Result) end proc; f_xx := proc (Functions::{Vector,
list({term, equation})}, Vars::{Vector, list}) local i1, i2, i3, f, x, fxx, nx
, nf, SetOfLHSs; if type(Functions,Vector) then f := convert(Functions,list) 
elif type(Functions,list(equation)) then f := map(rhs,Functions) else f := 
Functions end if; if type(Vars,Vector) then x := convert(Vars,list) else x :=
Vars end if; nx := nops(x); nf := nops(f); if not nops(f) = nops(x) then error
"dimensions of f and x are not the same" end if; fxx := array(1 .. nx,1 .. nx,
1 .. nx); for i1 to nx do for i2 to nx do for i3 to nx do fxx[i1,i2,i3] := 
diff(f[i1],x[i2],x[i3]) end do end do end do; return eval(fxx) end proc; f_xp
:= proc (Functions::{Vector, list({term, equation})}, Vars::{Vector, list}, 
Pars::{Vector, list}) local i1, i2, i3, f, x, p, fxp, nx, nf, np, SetOfLHSs; 
if type(Functions,Vector) then f := convert(Functions,list) elif type(
Functions,list(equation)) then f := map(rhs,Functions) else f := Functions end
if; if type(Vars,Vector) then x := convert(Vars,list) else x := Vars end if; 
if type(Pars,Vector) then p := convert(Pars,list) else p := Pars end if; nx :=
nops(x); nf := nops(f); np := nops(p); if not nops(f) = nops(x) then error 
"dimensions of f and x are not the same" end if; fxp := array(1 .. nx,1 .. nx,
1 .. np); for i1 to nx do for i2 to nx do for i3 to np do fxp[i1,i2,i3] := 
diff(f[i1],x[i2],p[i3]) end do end do end do; return eval(fxp) end proc; end 
module; FileOperations := module () export dirExists, fileExists, 
findFirstSharedObject, findObjectFile, nextFileNumber, removeAllFilesInDir, 
renameFile, runSystemCommandInDir, splitPathAndFile; dirExists := proc (
aString::string) local OldDir; OldDir := currentdir(); try currentdir(aString)
catch: RETURN(false) end try; currentdir(OldDir); RETURN(true) end proc; 
fileExists := proc (aString::string) local OldDir, IsDirect; OldDir := 
currentdir(); IsDirect := true; try currentdir(aString) catch: IsDirect := 
false end try; currentdir(OldDir); if FileTools[Exists](aString) and not 
IsDirect then RETURN(true) else RETURN(false) end if end proc; 
findFirstSharedObject := proc (TruncatedFileName::string) local MaxIndex, 
FileName, i1, SharedObjectEnd; MaxIndex := 999; SharedObjectEnd := ".so"; if 
_EnvOperatingSystemType = "windows" then SharedObjectEnd := ".dll" end if; for
i1 to MaxIndex do FileName := cat(TruncatedFileName,convert(i1,string),
SharedObjectEnd); if Aux:-FileOperations:-fileExists(FileName) then break end
if end do; if MaxIndex < i1 then error "requested file not found" end if; 
return FileName end proc; findObjectFile := proc (Filelist::list(string)) 
local OtherFiles, ObjectFiles, OldDir, item, Path, File; OtherFiles := []; 
ObjectFiles := []; OldDir := currentdir(); for item in Filelist while true do
Path, File := splitPathAndFile(item); try currentdir(Path); if fileExists(cat(
File,".o")) then ObjectFiles := [op(ObjectFiles), cat(Path,File,".o")] elif 
fileExists(cat(File,".obj")) then ObjectFiles := [op(ObjectFiles), cat(Path,
File,".obj")] else OtherFiles := [op(OtherFiles), item] end if catch: error 
finally currentdir(OldDir) end try end do; return ObjectFiles, OtherFiles end
proc; nextFileNumber := proc (aString::string) local i1, testForDir; i1 := 1;
 while fileExists(cat(aString,".",i1)) do i1 := i1+1 end do; RETURN(i1) end 
proc; removeAllFilesInDir := proc (aString::string) local SystemCommand, 
SystemCommandOk; SystemCommand := cat("rm -fr ",aString,"/* ",aString,
"/.[a-zA-Z0-9]*"); SystemCommandOk := system(SystemCommand); if not 
SystemCommandOk = 0 then error "removing all files in %1 was not successful",
aString end if end proc; renameFile := proc (OldFileName::string, NewFileName
::string) local UnixCommand; if not fileExists(OldFileName) then error 
"file to be renamed %1 does not exist", OldFileName end if; if fileExists(
NewFileName) then error "file named %1 already exists", NewFileName end if; 
UnixCommand := cat(`mv `,OldFileName,` `,NewFileName); system(UnixCommand) end
proc; runSystemCommandInDir := proc (command::string, WorkingDir::string) 
local OldDir; OldDir := currentdir(); try currentdir(WorkingDir); system(
command) catch: error finally currentdir(OldDir) end try; return true end proc
; splitPathAndFile := proc (PathAndFile::string) local SubsStrings, 
LeadingPath, FileName; if StringTools[Search]("/",PathAndFile) = 0 then return
"./", PathAndFile end if; SubsStrings := StringTools[Split](PathAndFile,"/");
LeadingPath := cat(seq(cat(SubsStrings[i1],"/"),i1 = 1 .. nops(SubsStrings)-1)
); FileName := SubsStrings[-1]; return LeadingPath, FileName end proc; end 
module; equalWithinRelErr := proc (x1::numeric, x2::numeric, relErr::numeric)
if abs(2*(x2-x1)/(x1+x2)) <= relErr then return true else return false end if
end proc; equalWithinFudgedRelErr := proc (x1::{numeric, list(numeric)}, x2::{
numeric, list(numeric)}, relErr::numeric) local i; if type(x1,numeric) and 
type(x2,numeric) then if getFudgedRelErr(x1,x2) <= relErr then return true 
else return false end if elif type(x1,list(numeric)) and type(x2,list(numeric)
) then if not nops(x1) = nops(x2) then error 
"lists specified in first and second parameter must have same length" end if;
for i to nops(x1) do if not equalWithinFudgedRelErr(x1[i],x2[i],relErr) then 
return false end if end do; return true else error 
"first and second parameter must be either type numeric or type list(numeric)"
end if end proc; getAbsErr := proc (x1::EvalsToFloat, x2::EvalsToFloat) option
inline; abs(x1-x2) end proc; getFudgedRelErr := proc (r1::numeric, r2::numeric
) local d, s, res; d := abs(r1-r2); s := abs(r1+r2); res := d/(1+s); return 
res end proc; getRelErr := proc (x1::EvalsToFloat, x2::EvalsToFloat) option 
inline; abs(2*evalf(x1-x2)/evalf(x1+x2)) end proc; init := proc () 
`type/DTASys` := proc (aDTASys) if Aux:-SystemClasses:-listOfErrorsInDTASys(
aDTASys) = [] then return true else return false end if end proc; 
`type/DAESys` := proc (aDAESys) if Aux:-SystemClasses:-listOfErrorsInDAESys(
aDAESys) = [] then return true else return false end if end proc; 
`type/EvalsToFloat` := proc (a::anything) return type(evalf(a),'float') end 
proc; `type/NLP` := proc (aSys) local ListOfErrors; ListOfErrors := Aux:-
SystemClasses:-listOfErrorsInNLP(aSys); if ListOfErrors = [] then return true
else return false end if end proc; `type/BCNLP` := proc (aNLP::table) if type(
eval(aNLP),NLP) and aNLP[Constraints] = [] and aNLP[LinearConstraints] = [] 
then return true else return false end if end proc; `type/term` := {`*`, `+`,
`^`, indexed, symbol, EvalsToFloat, function} end proc; IntervalArithmetics :=
module () export square; square := proc (Ix::[numeric, numeric]) local 
lowerBound, upperBound, IxBoundsSquared; IxBoundsSquared := Ix[1]^2, Ix[2]^2;
if Ix[1] <= 0 and 0 <= Ix[2] then lowerBound := 0 else lowerBound := min(
IxBoundsSquared) end if; upperBound := max(IxBoundsSquared); return [
lowerBound, upperBound] end proc; end module; LinearAlgebra := module () 
export calcGershgorinForRealSymMat; calcGershgorinForRealSymMat := proc (A::
Matrix) local numRows, numCols, D, i, j, lb, ub; numRows, numCols := :-
LinearAlgebra:-Dimension(A); if not numRows = numCols then error 
"first parameter must be square matrix" end if; D := Array(1 .. numRows); for
i to numRows do D[i] := 0; for j to numRows do if j = i then next end if; D[i]
:= D[i]+abs(A[i,j]) end do end do; lb := min(seq(A[i,i]-D[i],i = 1 .. numRows)
); ub := max(seq(A[i,i]+D[i],i = 1 .. numRows)); return [lb, ub] end proc; end
module; LinearEqns := module () local isLinearIn, convertLinSysToStandardForm;
export checkLinearityOfIn; isLinearIn := proc (ListOfRHS, Vars) local item, 
summand, ListOfIndets, SetOfIndets, VarsAsSet; VarsAsSet := convert(Vars,set);
for item in ListOfRHS while true do ListOfIndets := Aux:-ListOperations:-
getListOfIndetsIn(item); SetOfIndets := convert(ListOfIndets,set); if not 
`minus`(SetOfIndets,VarsAsSet) = {} then return false end if; if not (type(
item,linear(Vars)) or type(item,EvalsToFloat)) then return false end if end do
; return true end proc; convertLinSysToStandardForm := proc (ListOfRHS, Vars)
local i, LHS, A, b; A := linalg[jacobian](ListOfRHS,Vars); b := map(simplify,
evalm(evalm(`&*`(A,linalg[vector](Vars)))-linalg[vector](ListOfRHS))); return
evalm(A), evalm(b) end proc; checkLinearityOfIn := proc (Eq::{term, list(term)
, list({equation(term), term < term}), equation(term), term < term}, Vars::
list(name), A::name, b::name) local ListOfRHS; if type(Eq,list({equation, term
< term})) then ListOfRHS := [seq(rhs(Eq[i1])-lhs(Eq[i1]),i1 = 1 .. nops(Eq))]
elif type(Eq,list(term)) then ListOfRHS := Eq elif type(Eq,{equation(term), 
term < term}) then ListOfRHS := [rhs(Eq)-lhs(Eq)] else ListOfRHS := [Eq] end 
if; if isLinearIn(ListOfRHS,Vars) then if 3 < nargs then A, b := 
convertLinSysToStandardForm(ListOfRHS,Vars) elif 2 < nargs then WARNING(
"expecting two optional arguments") end if; return true else return false end
if end proc; end module; ListOperations := module () export 
getIndexToMinElement, getIndexToMaxElement, getIndicesToAllMaxElements, 
getIndicesToAllMinElements, getListOfIndetsIn, getListOfUndefExprIn, 
getMissingAndObsoleteNames, getNameFromDerivSymbol, getObsolExprInDAESys, 
getObsolExprInList, getPosOfLHSin, getPosOfObjectInList, getRHSofIn, 
getSetOfValidExprIn, isLHSin, prettyPrintListToFile, removeItemFromList, 
setRHSofInTo, subsEqnListIntoEqn, subsSubsList; getIndexToMinElement := proc (
l::list(numeric)) local i, indexToMin; indexToMin := 1; for i to nops(l) do if
l[i] < l[indexToMin] then indexToMin := i end if end do; return indexToMin end
proc; getIndexToMaxElement := proc (l::list(numeric)) local i, indexToMax; 
indexToMax := 1; for i to nops(l) do if l[indexToMax] < l[i] then indexToMax 
:= i end if end do; return indexToMax end proc; getIndicesToAllMaxElements :=
proc (l::list(numeric)) local maxValue, i, indicesToMaximumElements, newIndex;
maxValue := max(op(l)); indicesToMaximumElements := {}; for i to nops(l) do if
l[i] = maxValue then newIndex := i; indicesToMaximumElements := `union`(
indicesToMaximumElements,{newIndex}) end if end do; return 
indicesToMaximumElements end proc; getIndicesToAllMinElements := proc (l::list
(numeric)) local minValue, i, indicesToMinimumElements, newIndex; minValue :=
min(op(l)); indicesToMinimumElements := {}; for i to nops(l) do if l[i] = 
minValue then newIndex := i; indicesToMinimumElements := `union`(
indicesToMinimumElements,{newIndex}) end if end do; return 
indicesToMinimumElements end proc; getListOfIndetsIn := proc (anEqn::{anything
, equation}) local listOfIndets, i1, currentItemInList, newItemsOfList, 
ArgumentOfFunction, Exponent; if type(anEqn,equation) then listOfIndets := 
convert(indets(op(2,anEqn)),list) else listOfIndets := convert(indets(anEqn),
list) end if; i1 := 1;  while i1 <= nops(listOfIndets) do currentItemInList :=
op(i1,listOfIndets); newItemsOfList := []; if type(currentItemInList,function)
then ArgumentOfFunction := op(1,currentItemInList); newItemsOfList := convert(
indets(ArgumentOfFunction),list) elif type(currentItemInList,`^`) then 
ArgumentOfFunction := op(1,currentItemInList); Exponent := op(2,
currentItemInList); newItemsOfList := [op(convert(indets(ArgumentOfFunction),
list)), op(convert(indets(Exponent),list))] end if; if not newItemsOfList = []
then listOfIndets := subsop(i1 = NULL,listOfIndets); listOfIndets := [op(
listOfIndets), op(newItemsOfList)] else i1 := 1+i1 end if end do; listOfIndets
:= convert(convert(listOfIndets,set),list) end proc; getListOfUndefExprIn := 
proc (anEqn::{anything, equation}, aModelOrList::{list, list(equation)}) local
listOfIndets, AllowedExpressions, UndefExpr, i1; listOfIndets := 
getListOfIndetsIn(anEqn); AllowedExpressions := getSetOfValidExprIn(
aModelOrList); UndefExpr := []; for i1 to nops(listOfIndets) do if not member(
listOfIndets[i1],AllowedExpressions) then UndefExpr := [op(UndefExpr), 
listOfIndets[i1]] end if end do; RETURN(UndefExpr) end proc; 
getMissingAndObsoleteNames := proc (Names1::{list(name), list(string), set(
name), set(string)}, Names2::{list(name), list(string), set(name), set(string)
}) local SetNames1, SetNames2, Missing, Obsolete; SetNames1 := convert(Names1,
set); SetNames2 := convert(Names2,set); Missing := `minus`(SetNames2,SetNames1
); Obsolete := `minus`(SetNames1,SetNames2); return Missing, Obsolete end proc
; getNameFromDerivSymbol := proc (d::symbol) local dString, NewName; dString 
:= convert(d,string); if not dString[-1] = "'" then error 
"expecting first argument to end with character '" end if; NewName := convert(
dString[1 .. -2],name); return NewName end proc; getObsolExprInDAESys := proc
(aDAESys::table) local ListOfEquations, ListOfNames, ObsolExpr; ListOfNames :=
[op(aDAESys["DynVars"]), op(aDAESys["AlgVars"]), op(aDAESys["Parameters"]), op
(map(lhs,aDAESys["ExplicitAEs"]))]; ListOfEquations := [op(aDAESys[
"ExplicitAEs"]), op(aDAESys["AEs"]), op(aDAESys["ODEs"])]; ObsolExpr := 
getObsolExprInList(ListOfNames,ListOfEquations); RETURN(ObsolExpr) end proc; 
getObsolExprInList := proc (aListOfNames::{name, list({name, name = anything})
}, aListOfEquations::{equation, list(equation)}) local ListOfNames, 
ListOfEquations, ListOfRHSs, IndetsOfRHSs, ObsolExpr, i1, item; if type(
aListOfNames,name) then ListOfNames := [aListOfNames] elif type(aListOfNames,
list(name)) then ListOfNames := aListOfNames else ListOfNames := []; for item
in aListOfNames while true do if type(item,name) then ListOfNames := [op(
ListOfNames), item] else ListOfNames := [op(ListOfNames), lhs(item)] end if 
end do end if; if type(aListOfEquations,equation) then ListOfEquations := [
aListOfEquations] elif type(aListOfEquations,list(equation)) then 
ListOfEquations := aListOfEquations end if; ListOfRHSs := map(rhs,
ListOfEquations); IndetsOfRHSs := getListOfIndetsIn(ListOfRHSs); ObsolExpr :=
[]; for i1 to nops(ListOfNames) do if member(ListOfNames[i1],IndetsOfRHSs) 
then next else ObsolExpr := [op(ObsolExpr), ListOfNames[i1]] end if end do; 
RETURN(ObsolExpr) end proc; getPosOfLHSin := proc (aName::name, aList) local 
i1, PositionOfLHS, currentEqn; if not (2 < nargs and args[3] = ('NoChecks')) 
then if not isLHSin(aName,aList) then error `, 1st argument does not occur` 
end if; if not type(aList,{list({name, equation}), ('array')({name, equation})
}) then error "unexpected 2nd input argument type" end if end if; i1 := 0; 
PositionOfLHS := 0;  while PositionOfLHS = 0 and i1 < nops(aList) do i1 := i1+
1; currentEqn := aList[i1]; if type(currentEqn,equation) then if lhs(
currentEqn) = aName then PositionOfLHS := i1 end if else if currentEqn = aName
then PositionOfLHS := i1 end if end if end do; RETURN(PositionOfLHS) end proc;
getPosOfObjectInList := proc (item::{indexed, symbol, list({indexed, symbol})}
, l::list({indexed, symbol})) local i, numListItems, itemInString, 
currentItemInString, res; if type(item,list) then res := [seq(
getPosOfObjectInList(item[i],l),i = 1 .. nops(item))]; return res end if; 
numListItems := nops(l); itemInString := convert(item,string); for i to 
numListItems do currentItemInString := convert(l[i],string); if 
currentItemInString = itemInString then return i end if end do; error 
"unable to find requested item in list" end proc; getRHSofIn := proc (aname::
name, alist::{table, list(equation)}) local i1, rightHandSide, currentEqn; if
type(alist,table) then error 
`, do not use GetRHSofIn on tables, use indexed table instead` end if; if not
isLHSin(aname,alist) then error 
`, 1st argument does not occur in 2nd argument` end if; for i1 to nops(alist)
do currentEqn := op(i1,alist); if op(1,currentEqn) = aname then rightHandSide
:= op(2,currentEqn); break end if end do; return rightHandSide end proc; 
getSetOfValidExprIn := proc (aModelOrList::list({name, name = anything})) 
local ListOfValidExpr, SetOfValidExpr, item; ListOfValidExpr := []; for item 
in aModelOrList while true do if type(item,equation) then ListOfValidExpr := [
op(ListOfValidExpr), lhs(item)] else ListOfValidExpr := [op(ListOfValidExpr),
item] end if end do; SetOfValidExpr := convert(ListOfValidExpr,set); return 
SetOfValidExpr end proc; isLHSin := proc (aName::name, aList::{table, list({
name, equation}), ('array')({name, equation})}) local i1, eqn, NumEntries; if
aList = [] then RETURN(false) end if; if type(aList,{list({name, equation}), (
'array')({name, equation})}) then if type(aList,list) then NumEntries := nops(
aList) else NumEntries := linalg[vectdim](aList) end if; for i1 to NumEntries
do eqn := aList[i1]; if type(eqn,equation) then if lhs(eqn) = aName then 
RETURN(true) end if else if eqn = aName then RETURN(true) end if end if end do
elif type(aList,table) then if member([aName],[indices(aList)]) then RETURN(
true) end if end if; RETURN(false) end proc; prettyPrintListToFile := proc (
aList::list, listName::string, fileName::string) local numElements, i, line, 
fileHandle; numElements := nops(aList); if FileOperations:-fileExists(fileName
) then error cat("fileto be created ",fileName," already exists") else 
fileHandle := open(fileName,'WRITE') end if; if 3 < nargs then if type(args[4]
,list(string)) then for line in args[4] while true do fprintf(fileHandle,
"# %s\n",line) end do else error 
"expecting optional 4th parameter to be of type list(string)" end if end if; 
fprintf(fileHandle,"%s:= [\n",listName); for i to numElements-1 do fprintf(
fileHandle,"  %q,\n",aList[i]) end do; fprintf(fileHandle,"  %q\n",aList[
numElements]); fprintf(fileHandle,"];\n"); return end proc; removeItemFromList
:= proc (ReqItems::{name, string, list(name), list(string)}, aList::list({{
name, symbol} = anything, {string}, {name, symbol}})) local aName, ListOfNames
, item, item2, NewList, ItemsToBeRemoved; ListOfNames := []; for item in aList
while true do if type(item,equation) then ListOfNames := [op(ListOfNames), lhs
(item)] else ListOfNames := [op(ListOfNames), item] end if end do; if type(
ReqItems,name) then ItemsToBeRemoved := [ReqItems] else ItemsToBeRemoved := 
ReqItems end if; for aName in ItemsToBeRemoved while true do if not member(
aName,ListOfNames) then error 
"item %1 specified in first argument does not occur in second argument", aName
end if end do; NewList := []; for item2 in aList while true do if type(item2,
equation) then if not member(lhs(item2),ItemsToBeRemoved) then NewList := [op(
NewList), item2] end if else if not member(item2,ItemsToBeRemoved) then 
NewList := [op(NewList), item2] end if end if end do; return NewList end proc;
setRHSofInTo := proc (aName::name, FormerList::{table, list({name, equation})}
, NewRHS::anything) local NewList, PositionLHS; if not Aux:-ListOperations:-
isLHSin(aName,FormerList) then ERROR(
`, 1st argument does not occur in 2nd argument`) end if; if type(FormerList,
list({name, equation})) then PositionLHS := Aux:-ListOperations:-getPosOfLHSin
(aName,FormerList); NewList := FormerList; NewList := subsop(PositionLHS = (
aName = NewRHS),NewList) elif type(FormerList,table) then NewList := copy(
FormerList); NewList[aName] := NewRHS end if; return NewList end proc; 
subsEqnListIntoEqn := proc (aList::list(equation), anEqn::{list({`<`, term, 
equation}), {`<`, term, equation}}) local i1, i2, newEqn, NumberOfEqns, 
ListOfEqns, newList; NumberOfEqns := nops(aList); if type(anEqn,{`<`, term, 
equation}) then newEqn := subs(seq(aList[NumberOfEqns+1-i1],i1 = 1 .. 
NumberOfEqns),anEqn); RETURN(newEqn) elif type(anEqn,list({`<`, term, equation
})) then ListOfEqns := anEqn; newList := array(1 .. nops(ListOfEqns)); for i2
to nops(ListOfEqns) do newEqn := subs(seq(aList[NumberOfEqns+1-i1],i1 = 1 .. 
NumberOfEqns),ListOfEqns[i2]); newList[i2] := newEqn end do; newList := 
convert(newList,list); RETURN(eval(newList)) end if; false end proc; 
subsSubsList := proc (subsList::list(equation), target) local res, numEntries,
i; res := target; numEntries := nops(subsList); for i from numEntries by -1 to
1 do res := subs(subsList[i],res) end do; return res end proc; end module; NLP
:= module () export newNLP, parToVarInNLP, Scale, subsStandardNotationIntoNLP;
newNLP := proc () local aSys; aSys := table(); aSys["LinearConstraints"] := []
; aSys["Constraints"] := []; aSys["CostFunction"] := []; aSys["ExplicitAEs"] 
:= []; aSys["Parameters"] := []; aSys["Variables"] := []; return eval(aSys) 
end proc; parToVarInNLP := proc (ReqPars::list(name = EvalsToFloat .. 
EvalsToFloat), Sys::NLP) local NewNLP, item, NewVarNames, NewConstr, 
nNewConstr, NewLinConstr, nNewLinConstr, LinConstrForLinearityCheck, i1; 
NewNLP := copy(Sys); for item in ReqPars while true do NewNLP["Parameters"] :=
Aux:-ListOperations:-removeItemFromList(lhs(item),NewNLP["Parameters"]) end do
; NewNLP["Variables"] := [op(NewNLP["Variables"]), op(ReqPars)]; NewVarNames 
:= map(lhs,NewNLP["Variables"]); NewConstr := array(1 .. nops(Sys[
"LinearConstraints"])+nops(Sys["Constraints"])); nNewConstr := 0; NewLinConstr
:= array(1 .. nops(Sys["LinearConstraints"])); nNewLinConstr := 0; 
LinConstrForLinearityCheck := Aux:-ListOperations:-subsEqnListIntoEqn(Sys[
"ExplicitAEs"],Sys["LinearConstraints"]); for i1 to nops(Sys[
"LinearConstraints"]) do item := Sys["LinearConstraints"][i1]; if Aux:-
LinearEqns:-checkLinearityOfIn(LinConstrForLinearityCheck[i1],NewVarNames) 
then nNewLinConstr := nNewLinConstr+1; NewLinConstr[nNewLinConstr] := item 
else nNewConstr := nNewConstr+1; NewConstr[nNewConstr] := item end if end do;
NewNLP["LinearConstraints"] := convert(NewLinConstr,list); NewNLP[
"LinearConstraints"] := NewNLP["LinearConstraints"][1 .. nNewLinConstr]; 
NewConstr := convert(NewConstr,list); NewNLP["Constraints"] := [op(Sys[
"Constraints"]), op(NewConstr[1 .. nNewConstr])]; return eval(NewNLP) end proc
; Scale := module () export calcc, calcD, CreateInstanceForNLP, 
createSubsLists, replaceInfiniteRangesByTrivialRanges; calcc := proc (Ranges::
list(name = range)) local NumRanges, VarNames, a, b, c, i1; NumRanges := nops(
Ranges); VarNames := map(lhs,Ranges); c := array(1 .. NumRanges); for i1 to 
NumRanges do a := op(1,rhs(Ranges[i1])); b := op(2,rhs(Ranges[i1])); c[i1] :=
evalf((b+a)/(b-a)); if not type(c[i1],float) then error 
"entry %1 of vector c is not a float", c[i1] end if end do; return eval(c) end
proc; calcD := proc (Ranges::list(name = range)) local D, a, b, NumRanges, 
VarNames, i1; NumRanges := nops(Ranges); VarNames := map(lhs,Ranges); D := 
array(1 .. NumRanges); for i1 to NumRanges do a := op(1,rhs(Ranges[i1])); b :=
op(2,rhs(Ranges[i1])); D[i1] := evalf(2/(b-a)); if not type(D[i1],float) then
error "entry %1 on diagonal of diagonal matrix D is not a float", D[i1] end if
end do; return eval(D) end proc; CreateInstanceForNLP := proc (Sys::{list(name
= range), NLP}) local ReqRanges, Replaced, Ranges, VarNames, D, c, SubsLists,
SubsListPhysToScaled, SubsListScaledToPhys, InvD, i1, 
SubsListScaledParsPhysToScaled, SubsListScaledParsScaledToPhys, ScaledParNames
; if type(Sys,list) then ReqRanges := Sys else ReqRanges := Sys["Variables"] 
end if; Replaced, Ranges := replaceInfiniteRangesByTrivialRanges(ReqRanges); 
if not Replaced = [] then printf(
"replaced %d ranges which contained +/-infinity by trivial ranges -1..1\n",
nops(Replaced)) end if; VarNames := map(lhs,ReqRanges); D := calcD(Ranges); c
:= calcc(Ranges); SubsLists := createSubsLists(VarNames,D,c); 
SubsListPhysToScaled := SubsLists[1]; SubsListScaledToPhys := SubsLists[2]; 
InvD := array(1 .. nops(VarNames)); for i1 to nops(VarNames) do InvD[i1] := 1/
D[i1] end do; SubsListScaledParsPhysToScaled := []; 
SubsListScaledParsScaledToPhys := []; ScaledParNames := {}; module () local 
NamesScaledParsAndVars, NumVarRanges, NumScaledParRanges; export getc, getD, 
getInvD, getSubsListPhysToScaled, getSubsListScaledToPhys, getListOfVars, 
getMissingAssigns, getOutOfBounds, getUnknownAssigns, mapScaledToPhys, 
mapPhysToScaled, subsIntoDAESys, subsIntoExtAESys, subsIntoHopfDAESys, 
subsIntoNLP; NamesScaledParsAndVars := [op(VarNames), op(ScaledParNames)]; 
NumVarRanges := nops(VarNames); NumScaledParRanges := nops(ScaledParNames); 
getc := proc () return c end proc; getD := proc () return D end proc; getInvD
:= proc () return InvD end proc; getSubsListPhysToScaled := proc () local tmp;
tmp := [op(SubsListPhysToScaled), op(SubsListScaledParsPhysToScaled)]; return
tmp end proc; getSubsListScaledToPhys := proc () local tmp; tmp := [op(
SubsListScaledToPhys), op(SubsListScaledParsScaledToPhys)]; return tmp end 
proc; getListOfVars := proc () return VarNames end proc; getMissingAssigns :=
proc (SubsList::list(name = EvalsToFloat)) local Missing; Missing := `minus`(
convert(NamesScaledParsAndVars,set),convert(map(lhs,SubsList),set)); return 
Missing end proc; getOutOfBounds := proc (SubsList::list(name = EvalsToFloat))
local OutOfBounds, item; OutOfBounds := []; for item in SubsList while true do
if 1 < rhs(item) or rhs(item) < -1 then OutOfBounds := [op(OutOfBounds), lhs(
item)] end if end do; return OutOfBounds end proc; getUnknownAssigns := proc (
SubsList::list(name = EvalsToFloat)) local Unscaled, ValuesUnscaled; Unscaled
:= `minus`(convert(map(lhs,SubsList),set),convert(NamesScaledParsAndVars,set))
; ValuesUnscaled := [seq(Unscaled[i1] = subs(SubsList,Unscaled[i1]),i1 = 1 ..
nops(Unscaled))]; return ValuesUnscaled end proc; mapScaledToPhys := proc (
SubsList::list(name = EvalsToFloat)) local NewVars, i1, NamedListVars, x, p, 
Missing, OutOfBounds, ValuesUnknowns, result, NewPars, NamedListPars, xNamed,
pNamed; Missing := getMissingAssigns(SubsList); if not Missing = {} then error
"assignments are missing for %1", Missing end if; ValuesUnknowns := 
getUnknownAssigns(SubsList); x := subs(SubsList,VarNames); p := subs(SubsList,
ScaledParNames); xNamed := [seq(VarNames[i1] = x[i1],i1 = 1 .. nops(x))]; 
pNamed := [seq(ScaledParNames[i1] = p[i1],i1 = 1 .. nops(p))]; OutOfBounds :=
getOutOfBounds([op(xNamed), op(pNamed)]); if not OutOfBounds = [] then WARNING
("values for some variables and/or parameters are out of bounds [-1, 1]: %1",
OutOfBounds) end if; NewVars := array(1 .. NumVarRanges); for i1 to 
NumVarRanges do NewVars[i1] := InvD[i1]*(x[i1]+c[i1]) end do; NamedListVars :=
[seq(VarNames[i1] = NewVars[i1],i1 = 1 .. NumVarRanges)]; NewPars := array(1 
.. NumScaledParRanges); for i1 to NumScaledParRanges do NewPars[i1] := 
InvScaledParsD[i1]*(p[i1]+ScaledParsc[i1]) end do; NamedListPars := [seq(
ScaledParNames[i1] = NewPars[i1],i1 = 1 .. NumScaledParRanges)]; result := [op
(NamedListVars), op(NamedListPars), op(ValuesUnknowns)]; return result end 
proc; mapPhysToScaled := proc (SubsList::list(name = EvalsToFloat)) local 
NewVars, i1, NamedListVars, x, p, Missing, OutOfBounds, ValuesUnknowns, result
, NewPars, NamedListPars; Missing := getMissingAssigns(SubsList); if not 
Missing = {} then error "assignments are missing for %1", Missing end if; 
ValuesUnknowns := getUnknownAssigns(SubsList); x := subs(SubsList,VarNames); p
:= subs(SubsList,ScaledParNames); NewVars := array(1 .. NumVarRanges); for i1
to NumVarRanges do NewVars[i1] := D[i1]*x[i1]-c[i1] end do; NamedListVars := [
seq(VarNames[i1] = NewVars[i1],i1 = 1 .. NumVarRanges)]; NewPars := array(1 ..
NumScaledParRanges); for i1 to NumScaledParRanges do NewPars[i1] := 
ScaledParsD[i1]*p[i1]-ScaledParsc[i1] end do; NamedListPars := [seq(
ScaledParNames[i1] = NewPars[i1],i1 = 1 .. NumScaledParRanges)]; result := [op
(NamedListVars), op(NamedListPars), op(ValuesUnknowns)]; return result end 
proc; subsIntoDAESys := proc (aDAESys::DAESys) local NewSys, FormerODE, NewRhs
, i1; if not aDAESys["AEs"] = [] then error 
"use on algebraic equations currently not implemented" end if; if not aDAESys[
"DynVars"] = VarNames then warning("instance of NLP:-Scale can only be applied\
 to DAESys which has same \"DynVars\" as CreateInstances has received variable\
s") end if; NewSys := copy(aDAESys); NewSys["ExplicitAEs"] := array(1 .. nops(
aDAESys["ExplicitAEs"])); for i1 to nops(aDAESys["ExplicitAEs"]) do NewSys[
"ExplicitAEs"][i1] := subs(SubsListScaledToPhys,SubsListScaledParsScaledToPhys
,aDAESys["ExplicitAEs"][i1]) end do; NewSys["ExplicitAEs"] := convert(NewSys[
"ExplicitAEs"],list); NewSys["ODEs"] := array(1 .. nops(aDAESys["ODEs"])); for
i1 to nops(aDAESys["ODEs"]) do FormerODE := aDAESys["ODEs"][i1]; NewRhs := 
subs(SubsListScaledToPhys,SubsListScaledParsScaledToPhys,rhs(FormerODE)); 
NewRhs := D[i1]*NewRhs; NewSys["ODEs"][i1] := lhs(FormerODE) = NewRhs end do;
NewSys["ODEs"] := convert(NewSys["ODEs"],list); return eval(NewSys) end proc;
subsIntoExtAESys := proc (anExtAESys::ExtAESys) local NewSys, FormerODE, 
NewRhs, i1; if not anExtAESys["AEs"] = [] then error 
"use on algebraic equations currently not implemented" end if; if not 
anExtAESys["Variables"] = VarNames then error "instance of NLP:-Scale can only\
 be applied to \"ExtAESys\" which has same \"Variables\" as CreateInstances ha\
s received variables" end if; NewSys := copy(anExtAESys); NewSys["ExplicitAEs"
] := array(1 .. nops(anExtAESys["ExplicitAEs"])); for i1 to nops(anExtAESys[
"ExplicitAEs"]) do NewSys["ExplicitAEs"][i1] := subs(SubsListScaledToPhys,
SubsListScaledParsScaledToPhys,anExtAESys["ExplicitAEs"][i1]) end do; NewSys[
"ExplicitAEs"] := convert(NewSys["ExplicitAEs"],list); NewSys["Equations"] :=
array(1 .. nops(anExtAESys["Equations"])); for i1 to nops(anExtAESys[
"Equations"]) do FormerODE := anExtAESys["Equations"][i1]; NewRhs := subs(
SubsListScaledToPhys,SubsListScaledParsScaledToPhys,rhs(FormerODE)); NewRhs :=
D[i1]*NewRhs; NewSys["Equations"][i1] := lhs(FormerODE) = NewRhs end do; 
NewSys["Equations"] := convert(NewSys["Equations"],list); return eval(NewSys)
end proc; subsIntoNLP := proc (anNLP::NLP) local NewSys, NewRhs, i1, Ranges, 
Bounds, NumBounds, LowerBounds, UpperBounds, NumVars, NewBounds, CurrentName,
CurrentRange; if not map(lhs,anNLP["Variables"]) = VarNames then error "instan\
ce of NLP:-Scale can only be applied to NLP which has same \"Variables\" as Cr\
eateInstances has received variables" end if; NewSys := copy(anNLP); NewSys[
"ExplicitAEs"] := array(1 .. nops(anNLP["ExplicitAEs"])); for i1 to nops(anNLP
["ExplicitAEs"]) do NewSys["ExplicitAEs"][i1] := subs(SubsListScaledToPhys,
SubsListScaledParsScaledToPhys,anNLP["ExplicitAEs"][i1]) end do; NewSys[
"ExplicitAEs"] := convert(NewSys["ExplicitAEs"],list); NewSys[
"LinearConstraints"] := array(1 .. nops(anNLP["LinearConstraints"])); for i1 
to nops(anNLP["LinearConstraints"]) do NewSys["LinearConstraints"][i1] := subs
(SubsListScaledToPhys,SubsListScaledParsScaledToPhys,anNLP["LinearConstraints"
][i1]) end do; NewSys["LinearConstraints"] := convert(NewSys[
"LinearConstraints"],list); NewSys["Constraints"] := array(1 .. nops(anNLP[
"Constraints"])); for i1 to nops(anNLP["Constraints"]) do NewSys["Constraints"
][i1] := subs(SubsListScaledToPhys,SubsListScaledParsScaledToPhys,anNLP[
"Constraints"][i1]) end do; NewSys["Constraints"] := convert(NewSys[
"Constraints"],list); NewSys["CostFunction"] := subs(SubsListScaledToPhys,
SubsListScaledParsScaledToPhys,anNLP["CostFunction"]); NumVars := nops(anNLP[
"Variables"]); NewBounds := array(1 .. NumVars); for i1 to NumVars do 
CurrentName := lhs(anNLP["Variables"][i1]); CurrentRange := rhs(anNLP[
"Variables"][i1]); if op(1,CurrentRange) = -infinity or op(2,CurrentRange) = 
infinity then NewBounds[i1] := CurrentName = CurrentRange else NewBounds[i1] 
:= CurrentName = -1 .. 1 end if end do; NewSys["Variables"] := convert(
NewBounds,list); return eval(NewSys) end proc; subsIntoHopfDAESys := proc (
aDAESys::DAESys) local NewSys, FormerODE, NewRhs, i1; if not aDAESys["AEs"] =
[] then error "use on algebraic equations currently not implemented" end if; 
if not aDAESys["DynVars"] = VarNames[1 .. nops(aDAESys["DynVars"])] then 
warning("instance of LA:-Scale can only be applied to \"DAESys\" which has sam\
e \"DynVars\" as CreateInstances has received variables") end if; NewSys := 
copy(aDAESys); NewSys["ExplicitAEs"] := array(1 .. nops(aDAESys["ExplicitAEs"]
)); for i1 to nops(aDAESys["ExplicitAEs"]) do NewSys["ExplicitAEs"][i1] := 
subs(SubsListScaledToPhys,aDAESys["ExplicitAEs"][i1]) end do; NewSys[
"ExplicitAEs"] := convert(NewSys["ExplicitAEs"],list); NewSys["ODEs"] := array
(1 .. nops(aDAESys["ODEs"])); for i1 to nops(aDAESys["ODEs"]) do FormerODE :=
aDAESys["ODEs"][i1]; NewRhs := subs(SubsListScaledToPhys,rhs(FormerODE)); 
NewRhs := D[i1]*NewRhs; NewSys["ODEs"][i1] := lhs(FormerODE) = NewRhs end do;
NewSys["ODEs"] := convert(NewSys["ODEs"],list); return eval(NewSys) end proc 
end module end proc; createSubsLists := proc (VarNames, D::vector, c::vector)
local SubsListPhysToScaled, SubsListScaledToPhys, i1, NewItem, NumRanges; 
NumRanges := nops(VarNames); SubsListPhysToScaled := []; for i1 to NumRanges 
do NewItem := VarNames[i1] = D[i1]*VarNames[i1]-c[i1]; SubsListPhysToScaled :=
[op(SubsListPhysToScaled), NewItem] end do; SubsListScaledToPhys := []; for i1
to NumRanges do NewItem := VarNames[i1] = 1/D[i1]*(VarNames[i1]+c[i1]); 
SubsListScaledToPhys := [op(SubsListScaledToPhys), NewItem] end do; return [
SubsListPhysToScaled, SubsListScaledToPhys] end proc; 
replaceInfiniteRangesByTrivialRanges := proc (Ranges::list(name = range)) 
local NewRanges, item, LowerBound, UpperBound, ChangedRanges; NewRanges := 
Ranges; ChangedRanges := []; for item in NewRanges while true do LowerBound :=
op(1,rhs(item)); if LowerBound = -infinity then NewRanges := Aux:-
ListOperations:-setRHSofInTo(lhs(item),NewRanges,-1 .. 1); ChangedRanges := [
op(ChangedRanges), lhs(item)]; next end if; UpperBound := op(2,rhs(item)); if
UpperBound = infinity then NewRanges := Aux:-ListOperations:-setRHSofInTo(lhs(
item),NewRanges,-1 .. 1); ChangedRanges := [op(ChangedRanges), lhs(item)] end
if end do; return ChangedRanges, NewRanges end proc; end module; 
subsStandardNotationIntoNLP := proc (Sys::NLP, ReqNames::[name, name, name]) 
local NumOfVars, NumOfPars, ListOfVarsSubst, ListOfParsSubst, NewSystem, 
NamesVars, NamesPars, NumOfEAEs, NamesEAEs, ListOfEAEsSubst, StandardVarName,
StandardParName, StandardEAEName; if 1 < nargs then StandardParName := 
ReqNames[1]; StandardEAEName := ReqNames[2]; StandardVarName := ReqNames[3] 
else StandardParName := par; StandardEAEName := z; StandardVarName := x end if
; NumOfVars := nops(Sys["Variables"]); NumOfPars := nops(Sys["Parameters"]); 
NumOfEAEs := nops(Sys["ExplicitAEs"]); NamesVars := map(lhs,Sys["Variables"]);
ListOfVarsSubst := [seq(NamesVars[i1] = StandardVarName[i1],i1 = 1 .. 
NumOfVars)]; NamesPars := map(lhs,Sys["Parameters"]); ListOfParsSubst := [seq(
NamesPars[i1] = StandardParName[i1],i1 = 1 .. NumOfPars)]; NamesEAEs := map(
lhs,Sys["ExplicitAEs"]); ListOfEAEsSubst := [seq(NamesEAEs[i1] = 
StandardEAEName[i1],i1 = 1 .. NumOfEAEs)]; NewSystem := copy(Sys); NewSystem[
"Substitutions"] := [op(ListOfVarsSubst), op(ListOfParsSubst), op(
ListOfEAEsSubst)]; NewSystem["Constraints"] := subs(NewSystem["Substitutions"]
,Sys["Constraints"]); NewSystem["LinearConstraints"] := subs(NewSystem[
"Substitutions"],Sys["LinearConstraints"]); NewSystem["CostFunction"] := subs(
NewSystem["Substitutions"],Sys["CostFunction"]); NewSystem["ExplicitAEs"] := 
subs(NewSystem["Substitutions"],Sys["ExplicitAEs"]); NewSystem["Parameters"] 
:= subs(ListOfParsSubst,Sys["Parameters"]); NewSystem["Variables"] := subs(
ListOfVarsSubst,Sys["Variables"]); return eval(NewSystem) end proc; end module
; Other := module () export myModulo; myModulo := proc (n, m) local tmp; tmp 
:= `mod`(n,m); if tmp = 0 then return m else return tmp end if end proc; end 
module; Parsers := module () export readContentDataFromFile; 
readContentDataFromFile := proc (FileName::string) local line, ColNames, Data,
item, CurrentRow, i1, i2, alist, ScanString, ScanStringFirstLine, 
ListOfIndices, NumOfCols, NamedLists; try close(FileName) end try; line := 
readline(FileName); alist := convert(line,list); NumOfCols := 0; for i2 to 
nops(alist) do if alist[i2] = " " and not alist[i2+1] = " " then NumOfCols :=
NumOfCols+1 end if end do; if alist[1] = " " then NumOfCols = NumOfCols-1 end
if; NumOfCols := NumOfCols+1; ScanStringFirstLine := ""; ScanString := ""; for
i1 to NumOfCols do ScanStringFirstLine := cat(ScanStringFirstLine,"%a "); 
ScanString := cat(ScanString,"%e ") end do; ColNames := sscanf(line,
ScanStringFirstLine); ColNames := map(convert,ColNames,symbol); Data := table(
); for item in ColNames while true do Data[item] := [] end do;  do line := 
readline(FileName); if line = 0 then break end if; CurrentRow := sscanf(line,
ScanString); for i1 to nops(ColNames) do Data[ColNames[i1]] := [op(Data[
ColNames[i1]]), CurrentRow[i1]] end do end do; ListOfIndices := map(op,[
indices(Data)]); NamedLists := []; for i1 to nops(Data[ListOfIndices[1]]) do 
NamedLists := [op(NamedLists), [seq(ListOfIndices[i2] = Data[ListOfIndices[i2]
][i1],i2 = 1 .. nops(ListOfIndices))]] end do; return NamedLists end proc; end
module; Programming := module () export createHeader, printCompileOptions; 
createHeader := proc (procOrModuleName::string, author::string, fileName::
string) local copyrightIndicator, fh, today, copyrightNote; if 3 < nargs then
copyrightIndicator := args[4] else copyrightIndicator := "ISM" end if; if 
copyrightIndicator = "ISM" then copyrightNote := "# This software is developme\
nt software that is not free \n# and not open. It may become available in the \
future. Please\n# contact M. Moennigmann for details. \n# Copyright Integriert\
e Systeme der Mikroprozesstechnik." elif copyrightIndicator = "mmo" or 
copyrightIndicator = "Moennigmann" or copyrightIndicator = "Monnigmann" then 
copyrightNote := "# This software is development software that is not free \n#\
 and not open. It may become available in the future. Please\n# contact M. Moe\
nnigmann for details. Copyright M. Moennigmann." else error 
"internal error, unknown copyrightIndicator" end if; today := ssystem(
"date +%Y-%m-%d")[2][1 .. -2]; fh := open(fileName,'WRITE'); fprintf(fh,
"#------------------------------------------------------------\n"); fprintf(fh
,"#\n"); fprintf(fh,"# %s\n",procOrModuleName); fprintf(fh,"#\n"); fprintf(fh,
"# (put brief description here)\n"); fprintf(fh,"#\n"); fprintf(fh,"# Notes:\n\
"); fprintf(fh,"#\n"); fprintf(fh,"# to do:\n"); fprintf(fh,"#\n"); fprintf(fh
,"# revision history:\n"); fprintf(fh,"# %s Written by %s.\n",today,author); 
fprintf(fh,"#\n"); fprintf(fh,"%s\n",copyrightNote); fprintf(fh,"#\n"); 
fprintf(fh,"#------------------------------------------------------------\n");
fprintf(fh,"%s:= proc()\n",procOrModuleName); fprintf(fh,"end proc;\n"); close
(fh); return end proc; printCompileOptions := proc () local compileOptions; 
compileOptions := define_external('COMPILE_OPTIONS'); printf(
"compileOptions:-COMPILER:=          %q\n",compileOptions:-COMPILER); printf(
"compileOptions:-CFLAGS:=            %q\n",compileOptions:-CFLAGS); printf(
"compileOptions:-COMPILE_ONLY_FLAG:= %q\n",compileOptions:-COMPILE_ONLY_FLAG);
printf("compileOptions:-INC_FLAG:=          %q\n",compileOptions:-INC_FLAG); 
printf("compileOptions:-COBJ_FLAG:=         %q\n",compileOptions:-COBJ_FLAG);
printf("compileOptions:-FILE:=              %q\n",compileOptions:-FILE); 
printf("compileOptions:-FILE_EXT:=          %q\n",compileOptions:-FILE_EXT); 
printf("compileOptions:-INC_PATH:=          %q\n",compileOptions:-INC_PATH); 
printf("compileOptions:-OBJ_EXT:=           %q\n",compileOptions:-OBJ_EXT); 
printf("compileOptions:-COMPILE_COMMAND:=   %q\n",compileOptions:-
COMPILE_COMMAND); printf("compileOptions:-LINKER:=            %q\n",
compileOptions:-LINKER); printf("compileOptions:-LINK_FLAGS:=        %q\n",
compileOptions:-LINK_FLAGS); printf("compileOptions:-LIB_PATH:=          %q\n"
,compileOptions:-LIB_PATH); printf("compileOptions:-DLL_EXT:=           %q\n",
compileOptions:-DLL_EXT); printf("compileOptions:-SYS_LIBS:=          %q\n",
compileOptions:-SYS_LIBS); printf("compileOptions:-LIB:=               %q\n",
compileOptions:-LIB); printf("compileOptions:-LIBS:=              %q\n",
compileOptions:-LIBS); printf("compileOptions:-LIB_FLAG:=          %q\n",
compileOptions:-LIB_FLAG); printf("compileOptions:-LOBJ_FLAG:=         %q\n",
compileOptions:-LOBJ_FLAG); printf("compileOptions:-EXPORT_FLAG:=       %q\n",
compileOptions:-EXPORT_FLAG); printf("compileOptions:-FUNCTION:=          %q\n\
",compileOptions:-FUNCTION); printf("compileOptions:-LINK_COMMAND:=      %q\n"
,compileOptions:-LINK_COMMAND); return end proc; end module; SystemClasses :=
module () export exportDTASysToCLMatContM, exportDTASysToMatCont, 
listOfErrorsInDAESys, listOfErrorsInDAESysPart2, listOfErrorsInDTASys, 
listOfErrorsInDTASysPart2, listOfErrorsInEAEs, listOfErrorsInNlpPart2, 
listOfErrorsInNLP, noNameConflictsInUnitsInDAESys, 
noNameConflictsInUnitsInDTASys, subsExplicitAEsIntoDAESys, 
subsExplicitAlgEqnsIntoDTASys; exportDTASysToCLMatContM := proc (nameOfSys::
string, aSys::table) local DESys, headOfFile, paramOfSys, strParam, setForSubs
, dynVars, i, i1, i2, dynEqns, strDynEqns, funcDyDtOfFile, y0, strY0, 
initOfFile, lastPartOfFile, finalText; DESys := Aux:-SystemClasses:-
subsExplicitAlgEqnsIntoDTASys(aSys); headOfFile := StringTools[Join]([
"function out = ", convert(nameOfSys,string), "\nout{1} = @init;\nout{2} = @fu\
n_eval;\nout{3} = [];\nout{4} = [];\nout{5} = [];\nout{6} = [];\nout{7} = [];\
\nout{8} = [];\nout{9} = [];"],""); paramOfSys := map(lhs,DESys["Parameters"])
; strParam := ""; for i to nops(DESys["Parameters"]) do strParam := 
StringTools[Join]([strParam, convert(paramOfSys[i],string)],",") end do; 
setForSubs := {}; dynVars := DESys["DynVars"]; for i1 to nops(dynVars) do 
setForSubs := `union`(setForSubs,{dynVars[i1] = kmrgd(i1)}) end do; dynEqns :=
subs(setForSubs,map(rhs,DESys["DynEqns"])); strDynEqns := ""; for i2 to nops(
dynEqns) do strDynEqns := StringTools[Join]([convert(dynEqns[i2],string), 
strDynEqns],";;") end do; funcDyDtOfFile := StringTools[Join]([
"function dydt = fun_eval(t,kmrgd", strParam, ")\n", "dydt =[", strDynEqns, 
"];"],""); y0 := [seq(0,i = 1 .. nops(DESys["DynVars"]))]; strY0 := convert(y0
,string); initOfFile := StringTools[Join]([
"function [tspan,y0,options] = init\nhandles = feval(", nameOfSys, ");\ny0=",
strY0, ";\noptions = odeset('Jacobian',[],'JacobianP',[],'Hessians',[],'Hessia\
nsP',[]);\ntspan = [0 10];"],""); lastPartOfFile := StringTools[Join]([
"function jac = jacobian(t,kmrgd", strParam, 
")\nfunction jacp = jacobianp(t,kmrgd", strParam, 
")\nfunction hess = hessians(t,kmrgd", strParam, 
")\nfunction hessp = hessiansp(t,kmrgd", strParam, 
")\nfunction tens3  = der3(t,kmrgd ", strParam, 
")\nfunction tens4  = der4(t,kmrgd ", strParam, 
")\n function tens5  = der5(t,kmrgd", strParam, ")"],""); finalText := 
StringTools[Join]([headOfFile, "\n\n", funcDyDtOfFile, "\n\n", initOfFile, 
"\n\n", lastPartOfFile],""); RETURN(finalText) end proc; exportDTASysToMatCont
:= proc (aSys::table) local EqnsToBeExported, DESys, eqns, i; DESys := Aux:-
SystemClasses:-subsExplicitAlgEqnsIntoDTASys(aSys); eqns := DESys["DynEqns"];
EqnsToBeExported := eqns; for i to nops(eqns) do EqnsToBeExported[i] := 
convert(StringTools[Join]([convert(lhs(eqns[i]),string), "'"],""),name) = rhs(
eqns[i]) end do; RETURN([codegen[C](EqnsToBeExported), cat("Coordinates=",
DESys["DynVars"]), cat("Parameters=",map(lhs,DESys["Parameters"]))]) end proc;
listOfErrorsInDAESys := proc (aDAESys) local i1, i, ParsAndIndependents, 
SetOfValidExpr, NumberOfExplicitAEs, ListOfErrors, SetOfIndices, 
ListOfUndefExpr, LHSidesOfODEs, ObsolExpr, ParAndVar, DerivativesOfDynVars, 
ParAndExplAE, Integrators, item; if not type(aDAESys,table) then return false
end if; ListOfErrors := []; SetOfIndices := {indices(aDAESys)}; if not member(
["ExplicitAEs"],SetOfIndices) then ListOfErrors := [op(ListOfErrors), 
"entry ExplicitAEs missing"]; RETURN(ListOfErrors) elif not type(aDAESys[
"ExplicitAEs"],list(equation)) then ListOfErrors := [op(ListOfErrors), 
"entry ExplicitAEs must be a list of equations"]; RETURN(ListOfErrors) end if;
if not member(["ODEs"],SetOfIndices) then ListOfErrors := [op(ListOfErrors), 
"entry ODEs missing"]; RETURN(ListOfErrors) elif not type(aDAESys["ODEs"],list
(equation)) then ListOfErrors := [op(ListOfErrors), 
"entry ODEs must be a list of equations"]; RETURN(ListOfErrors) end if; if not
member(["Parameters"],SetOfIndices) then ListOfErrors := [op(ListOfErrors), 
"entry Parameters missing"]; RETURN(ListOfErrors) elif not type(aDAESys[
"Parameters"],list(name = EvalsToFloat)) then ListOfErrors := [op(ListOfErrors
), "entry Parameters must be a list of name= EvalsToFloat"]; RETURN(
ListOfErrors) end if; if not member(["DynVars"],SetOfIndices) then 
ListOfErrors := [op(ListOfErrors), "entry DynVars missing"]; RETURN(
ListOfErrors) elif not type(aDAESys["DynVars"],list) then ListOfErrors := [op(
ListOfErrors), "entry DynVars must be a list"]; RETURN(ListOfErrors) end if; 
if not member(["AEs"],SetOfIndices) then ListOfErrors := [op(ListOfErrors), 
"entry AEs missing"]; RETURN(ListOfErrors) elif not type(aDAESys["AEs"],list)
then ListOfErrors := [op(ListOfErrors), "entry AEs must be a list"]; RETURN(
ListOfErrors) end if; if not member(["AlgVars"],SetOfIndices) then 
ListOfErrors := [op(ListOfErrors), "entry AlgVars missing"]; RETURN(
ListOfErrors) elif not type(aDAESys["AlgVars"],list) then ListOfErrors := [op(
ListOfErrors), "entry AlgVars must be a list"]; RETURN(ListOfErrors) end if; 
if not ListOfErrors = [] then return ListOfErrors end if; if 1 < nargs then if
args[2] = ('long') or args[2] = ('all') or args[2] = ('strict') then 
ListOfErrors := listOfErrorsInDAESysPart2(aDAESys) end if end if; return 
ListOfErrors end proc; listOfErrorsInDAESysPart2 := proc (aDAESys) local 
ListOfErrors, LHSidesOfODEs, DynVarsFromODEs, i, ParsAndIndependents, 
NumberOfExplicitAEs, SetOfValidExpr, i1, ListOfUndefExpr, ObsolExpr, 
Integrators, item, ParAndVar, ParAndExplAE; ListOfErrors := []; LHSidesOfODEs
:= map(lhs,aDAESys["ODEs"]); DynVarsFromODEs := map(Aux:-ListOperations:-
getNameFromDerivSymbol,LHSidesOfODEs); DynVarsFromODEs := map(convert,
DynVarsFromODEs,string); if not DynVarsFromODEs = map(convert,aDAESys[
"DynVars"],string) then ListOfErrors := [op(ListOfErrors), 
" DynVars must be in same order as lhs of ODEs"] end if; if not nops(aDAESys[
"AEs"]) = nops(aDAESys["AlgVars"]) then ListOfErrors := [op(ListOfErrors), 
`number of algebraic variables is not equal to number of AEs`] end if; for i 
to nops(aDAESys["AEs"]) do if not lhs(aDAESys["AEs"][i]) = 0 then ListOfErrors
:= [op(ListOfErrors), `left hand sides of AEs must be ZERO`] end if end do; 
ParsAndIndependents := `union`(`union`(Aux:-ListOperations:-
getSetOfValidExprIn(aDAESys["Parameters"]),Aux:-ListOperations:-
getSetOfValidExprIn(aDAESys["AlgVars"])),Aux:-ListOperations:-
getSetOfValidExprIn(aDAESys["DynVars"])); NumberOfExplicitAEs := nops(aDAESys[
"ExplicitAEs"]); SetOfValidExpr := ParsAndIndependents; for i1 to 
NumberOfExplicitAEs do ListOfUndefExpr := Aux:-ListOperations:-
getListOfUndefExprIn(aDAESys["ExplicitAEs"][i1],convert(SetOfValidExpr,list));
if not ListOfUndefExpr = [] then ListOfErrors := [op(ListOfErrors), cat(
ListOfUndefExpr,` in equation `,aDAESys["ExplicitAEs"][i1],
` of ExplicitAEs undefined.`)] end if; SetOfValidExpr := `union`(
SetOfValidExpr,{op(1,aDAESys["ExplicitAEs"][i1])}) end do; for i1 to nops(
aDAESys["ODEs"]) do ListOfUndefExpr := Aux:-ListOperations:-
getListOfUndefExprIn(aDAESys["ODEs"][i1],convert(SetOfValidExpr,list)); if not
ListOfUndefExpr = [] then ListOfErrors := [op(ListOfErrors), cat(
ListOfUndefExpr,` in equation `,aDAESys["ODEs"][i1],` of ODEs undefined.`)] 
end if end do; for i1 to nops(aDAESys["AEs"]) do ListOfUndefExpr := Aux:-
ListOperations:-getListOfUndefExprIn(aDAESys["AEs"][i1],convert(SetOfValidExpr
,list)); if not ListOfUndefExpr = [] then ListOfErrors := [op(ListOfErrors), 
cat(ListOfUndefExpr,` in equation `,aDAESys["AEs"][i1],` of AEs undefined.`)]
end if end do; ObsolExpr := Aux:-ListOperations:-getObsolExprInDAESys(aDAESys)
; Integrators := []; for item in ObsolExpr while true do if member(item,
aDAESys["DynVars"]) then Integrators := [op(Integrators), item] end if end do;
if not Integrators = [] then ObsolExpr := Aux:-ListOperations:-
removeItemFromList(Integrators,ObsolExpr); WARNING("in ListOfErrorsInDAESys, t\
he following variables only occur as derivatives but not as variables: %1",op(
Integrators)) end if; if not ObsolExpr = [] then WARNING(
"in ListOfErrorsInDAESys, %1 is an obsolete expression",op(ObsolExpr)) end if;
ParAndVar := map(has,aDAESys["Parameters"],aDAESys["DynVars"]); if member(true
,ParAndVar,'i2') then ListOfErrors := [op(ListOfErrors), cat(aDAESys[
"Parameters"][i2],` is declared both parameter and variable `)] end if; 
ParAndExplAE := map(has,aDAESys["Parameters"],map(lhs,aDAESys["ExplicitAEs"]))
; if member(true,ParAndExplAE,'i2') then ListOfErrors := [op(ListOfErrors), 
cat(aDAESys["Parameters"][i2],
` occurs both in Parameters and lhs of ExplicitAEs `)] end if; if member('
Units',map(op,[indices(aDAESys)])) then noNameConflictsInUnitsInDAESys(aDAESys
) end if; return ListOfErrors end proc; listOfErrorsInDTASys := proc (aDTASys)
local i1, i, ParsAndIndependents, SetOfValidExpr, NumberOfExplicitAlgEqns, 
ListOfErrors, SetOfIndices, ListOfUndefExpr, LHSidesOfODEs, ObsolExpr, 
ParAndVar, DerivativesOfDynVars, ParAndExplAE, Integrators, item; if not type(
aDTASys,table) then return false end if; ListOfErrors := []; SetOfIndices := {
indices(aDTASys)}; if not member(["ExplicitAlgEqns"],SetOfIndices) then 
ListOfErrors := [op(ListOfErrors), "entry \"ExplicitAlgEqns\" missing"]; 
RETURN(ListOfErrors) elif not type(aDTASys["ExplicitAlgEqns"],list(equation))
then ListOfErrors := [op(ListOfErrors), 
"entry \"ExplicitAlgEqns\" must be a list of equations"]; RETURN(ListOfErrors)
end if; if not member(["DynEqns"],SetOfIndices) then ListOfErrors := [op(
ListOfErrors), "entry \"DynEqns\" missing"]; RETURN(ListOfErrors) elif not 
type(aDTASys["DynEqns"],list(equation)) then ListOfErrors := [op(ListOfErrors)
, "entry \"DynEqns\" must be a list of equations"]; RETURN(ListOfErrors) end 
if; if not member(["Parameters"],SetOfIndices) then ListOfErrors := [op(
ListOfErrors), "entry \"Parameters\" missing"]; RETURN(ListOfErrors) elif not
type(aDTASys["Parameters"],list(name = EvalsToFloat)) then ListOfErrors := [op
(ListOfErrors), "entry \"Parameters\" must be a list of name= EvalsToFloat"];
RETURN(ListOfErrors) end if; if not member(["DynVars"],SetOfIndices) then 
ListOfErrors := [op(ListOfErrors), "entry \"DynVars\" missing"]; RETURN(
ListOfErrors) elif not type(aDTASys["DynVars"],list) then ListOfErrors := [op(
ListOfErrors), "entry \"DynVars\" must be a list"]; RETURN(ListOfErrors) end 
if; if not member(["AlgEqns"],SetOfIndices) then ListOfErrors := [op(
ListOfErrors), "entry \"AlgEqns\" missing"]; RETURN(ListOfErrors) elif not 
type(aDTASys["AlgEqns"],list) then ListOfErrors := [op(ListOfErrors), 
"entry \"AlgEqns\" must be a list"]; RETURN(ListOfErrors) end if; if not 
member(["AlgVars"],SetOfIndices) then ListOfErrors := [op(ListOfErrors), 
"entry \"AlgVars\" missing"]; RETURN(ListOfErrors) elif not type(aDTASys[
"AlgVars"],list) then ListOfErrors := [op(ListOfErrors), 
"entry \"AlgVars\" must be a list"]; RETURN(ListOfErrors) end if; if not 
ListOfErrors = [] then return ListOfErrors end if; if 1 < nargs then if args[2
] = ('long') or args[2] = ('all') or args[2] = ('strict') then ListOfErrors :=
listOfErrorsInDTASysPart2(aDTASys) end if end if; return ListOfErrors end proc
; listOfErrorsInDTASysPart2 := proc (aDTASys) local ListOfErrors, 
LHSidesOfDynEqns, DynVarsFromDynEqns, i, ParsAndIndependents, 
NumberOfExplicitAlgEqns, SetOfValidExpr, i1, ListOfUndefExpr, ObsolExpr, 
Integrators, item, ParAndVar, ParAndExplAE, DynVarsAndExplAE, aListOfNames, 
aListOfEquations; ListOfErrors := []; DynVarsFromDynEqns := map(lhs,aDTASys[
"DynEqns"]); DynVarsFromDynEqns := map(convert,DynVarsFromDynEqns,string); if
not DynVarsFromDynEqns = map(convert,aDTASys["DynVars"],string) then 
ListOfErrors := [op(ListOfErrors), 
" \"DynVars\" must be in same order as lhs of \"DynEqns\" "] end if; if not 
nops(aDTASys["AlgEqns"]) = nops(aDTASys["AlgVars"]) then ListOfErrors := [op(
ListOfErrors), 
`number of algebraic variables is not equal to number of AlgEqns`] end if; for
i to nops(aDTASys["AlgEqns"]) do if not lhs(aDTASys["AlgEqns"][i]) = 0 then 
ListOfErrors := [op(ListOfErrors), `left hand sides of "AlgEqns" must be ZERO`
] end if end do; ParAndVar := map(has,aDTASys["Parameters"],aDTASys["DynVars"]
); if member(true,ParAndVar,'i2') then ListOfErrors := [op(ListOfErrors), cat(
aDTASys["Parameters"][i2],` is declared both in "Parameters" and "DynVars" `)]
end if; ParAndExplAE := map(has,aDTASys["Parameters"],map(lhs,aDTASys[
"ExplicitAlgEqns"])); if member(true,ParAndExplAE,'i2') then ListOfErrors := [
op(ListOfErrors), cat(aDTASys["Parameters"][i2],
` occurs both in "Parameters" and lhs of "ExplicitAlgEqns" `)] end if; 
DynVarsAndExplAE := map(has,aDTASys["DynVars"],map(lhs,aDTASys[
"ExplicitAlgEqns"])); if member(true,DynVarsAndExplAE,'i2') then ListOfErrors
:= [op(ListOfErrors), cat(aDTASys["DynVars"][i2],
` occurs both in "DynVars" and lhs of "ExplicitAlgEqns" `)] end if; 
ParsAndIndependents := `union`(`union`(Aux:-ListOperations:-
getSetOfValidExprIn(aDTASys["Parameters"]),Aux:-ListOperations:-
getSetOfValidExprIn(aDTASys["AlgVars"])),Aux:-ListOperations:-
getSetOfValidExprIn(aDTASys["DynVars"])); NumberOfExplicitAlgEqns := nops(
aDTASys["ExplicitAlgEqns"]); SetOfValidExpr := ParsAndIndependents; for i1 to
NumberOfExplicitAlgEqns do ListOfUndefExpr := Aux:-ListOperations:-
getListOfUndefExprIn(aDTASys["ExplicitAlgEqns"][i1],convert(SetOfValidExpr,
list)); if not ListOfUndefExpr = [] then ListOfErrors := [op(ListOfErrors), 
cat(ListOfUndefExpr,` in equation `,aDTASys["ExplicitAlgEqns"][i1],
` of "ExplicitAlgEqns" undefined.`)] end if; SetOfValidExpr := `union`(
SetOfValidExpr,{op(1,aDTASys["ExplicitAlgEqns"][i1])}) end do; for i1 to nops(
aDTASys["DynEqns"]) do ListOfUndefExpr := Aux:-ListOperations:-
getListOfUndefExprIn(aDTASys["DynEqns"][i1],convert(SetOfValidExpr,list)); if
not ListOfUndefExpr = [] then ListOfErrors := [op(ListOfErrors), cat(
ListOfUndefExpr,` in equation `,aDTASys["DynEqns"][i1],
` of "DynEqns" undefined.`)] end if end do; aListOfNames := [op(aDTASys[
"AlgVars"]), op(aDTASys["Parameters"]), op(map(lhs,aDTASys["ExplicitAlgEqns"])
)]; aListOfEquations := [op(aDTASys["ExplicitAlgEqns"]), op(aDTASys["AlgEqns"]
), op(aDTASys["DynEqns"])]; ObsolExpr := Aux:-ListOperations:-
getObsolExprInList(aListOfNames,aListOfEquations); if not ObsolExpr = [] then
ListOfErrors := [op(ListOfErrors), cat(ObsolExpr,` is defined but not used.`)]
end if; if not nops({op(map(lhs,aDTASys["Parameters"]))}) = nops([op(map(lhs,
aDTASys["Parameters"]))]) then ListOfErrors := [op(ListOfErrors), cat([op(map(
lhs,aDTASys["Parameters"]))],`some of "Parameters" is defined multiple times.`
)] end if; if not nops({op(map(lhs,aDTASys["ExplicitAlgEqns"]))}) = nops([op(
map(lhs,aDTASys["ExplicitAlgEqns"]))]) then ListOfErrors := [op(ListOfErrors),
cat([op(map(lhs,aDTASys["ExplicitAlgEqns"]))],
`some of lhs of "ExplicitAlgEqns" is defined multiple times.`)] end if; if 
member('Units',map(op,[indices(aDTASys)])) then noNameConflictsInUnitsInDTASys
(aDTASys) end if; return ListOfErrors end proc; listOfErrorsInEAEs := proc (
EAEs::list(name = term), ReqValidSymbs::{list(name), set(name)}) local 
ValidSymbs, NumEAEs, FoundNoErrors, i1, Inds, i2; ValidSymbs := convert(
ReqValidSymbs,set); NumEAEs := nops(EAEs); FoundNoErrors := true; for i1 to 
NumEAEs do if `mod`(i1,100) = 0 then printf(
"have analysed %d out of %d equations\n",i1,NumEAEs) end if; Inds := convert(
Aux:-ListOperations:-getListOfIndetsIn(rhs(EAEs[i1])),set); Inds := `minus`(
Inds,ValidSymbs); for i2 from i1-1 by -1 to 1 do Inds := `minus`(Inds,{lhs(
EAEs[i2])}); if Inds = {} then break end if end do; if not Inds = {} then 
FoundNoErrors := false; printf("symbols %q \n",Inds); printf(
"   cannot be resolved in equation number %d\n",i1) end if end do; return 
FoundNoErrors end proc; listOfErrorsInNlpPart2 := proc (Sys) local 
ListOfErrors, i1, NumberOfConstraints, NumberOfLinearConstraints, 
SetOfValidExpr, ListOfValidExpr, ListOfUndefExpr, ParAndVar, ObsolExpr, 
ListOfTerms, item, SetOfVariables, SetOfParameters, ListOfRhs, Rhs, 
IndetsOfConstraints; SetOfVariables := Aux:-ListOperations:-
getSetOfValidExprIn(map(lhs,Sys["Variables"])); SetOfParameters := Aux:-
ListOperations:-getSetOfValidExprIn(map(lhs,Sys["Parameters"])); 
SetOfValidExpr := `union`(SetOfParameters,SetOfVariables); ListOfValidExpr :=
convert(SetOfValidExpr,list); ListOfErrors := listOfErrorsInEAEs(Sys[
"ExplicitAEs"],ListOfValidExpr); ListOfValidExpr := [op(ListOfValidExpr), op(
map(lhs,Sys["ExplicitAEs"]))]; NumberOfConstraints := nops(Sys["Constraints"])
; for i1 to NumberOfConstraints do ListOfUndefExpr := Aux:-ListOperations:-
getListOfUndefExprIn(Sys["Constraints"][i1],ListOfValidExpr); if not 
ListOfUndefExpr = [] then ListOfErrors := [op(ListOfErrors), cat(
"expression(s) ",ListOfUndefExpr," undefined in constraint number ",i1,
" of NLP[\"Constraints\"]")] end if end do; NumberOfLinearConstraints := nops(
Sys["LinearConstraints"]); for i1 to NumberOfLinearConstraints do 
ListOfUndefExpr := Aux:-ListOperations:-getListOfUndefExprIn(Sys[
"LinearConstraints"][i1],ListOfValidExpr); if not ListOfUndefExpr = [] then 
ListOfErrors := [op(ListOfErrors), cat("expression(s) ",ListOfUndefExpr,
" undefined in constraint number ",i1," of NLP[\"LinearConstraints\"]")] end 
if end do; ListOfUndefExpr := Aux:-ListOperations:-getListOfUndefExprIn(Sys[
"CostFunction"],ListOfValidExpr); if not ListOfUndefExpr = [] then 
ListOfErrors := [op(ListOfErrors), cat("expression(s) ",ListOfUndefExpr,
" undefined in CostFunction")] end if; ParAndVar := [op(`intersect`(
SetOfParameters,SetOfVariables))]; if not ParAndVar = [] then ListOfErrors :=
[op(ListOfErrors), cat(ParAndVar,
` is/are declared both, parameter and variable `)] end if; if not Sys[
"LinearConstraints"] = [] then ListOfRhs := Aux:-ListOperations:-
subsEqnListIntoEqn(Sys["ExplicitAEs"],Sys["LinearConstraints"]); ListOfRhs :=
map(rhs,ListOfRhs); if not Aux:-LinearEqns:-checkLinearityOfIn(ListOfRhs,map(
lhs,Sys["Variables"])) then ListOfErrors := [op(ListOfErrors), 
"Not all constraints in LinearConstraints are linear"] end if end if; return 
ListOfErrors end proc; listOfErrorsInNLP := proc (Sys) local ListOfErrors, 
SetOfMissingIndices; if not type(Sys,table) then ListOfErrors := [
"instance of NLP is expected to be of type table"]; return ListOfErrors end if
; ListOfErrors := []; SetOfMissingIndices := `minus`({["Constraints"], [
"CostFunction"], ["ExplicitAEs"], ["Parameters"], ["Variables"], [
"LinearConstraints"]},{indices(Sys)}); if not SetOfMissingIndices = {} then 
ListOfErrors := [op(ListOfErrors), cat("entries",op(SetOfMissingIndices),
"missing in NLP")]; RETURN(ListOfErrors) end if; if Sys["LinearConstraints"] =
[] and Sys["Constraints"] = [] and Sys["CostFunction"] = [] then ListOfErrors
:= [op(ListOfErrors), cat(
"at least one \"LinearConstraint\" or one \"Constraint\" or the ",
"\"CostFunction\" must be given")]; RETURN(ListOfErrors) end if; if Sys[
"Variables"] = [] then ListOfErrors := [op(ListOfErrors), 
"List of variables must not be empty"]; RETURN(ListOfErrors) end if; if not 
type(Sys["ExplicitAEs"],list(name = term)) then ListOfErrors := [op(
ListOfErrors), cat("entry \"ExplicitAEs\" must be a list with expressions ",
"of form 0<term or 0=term; it may be empty")]; RETURN(ListOfErrors) end if; if
not type(Sys["LinearConstraints"],{list({0 = term, 0 < term}), []}) then 
ListOfErrors := [op(ListOfErrors), cat(
"entry \"LinearConstraints\" must be a list with expressions ",
"of form 0<term or 0=term; it may be empty")]; RETURN(ListOfErrors) end if; if
not type(Sys["Constraints"],list({0 = term, 0 < term})) then ListOfErrors := [
op(ListOfErrors), cat(
"entry \"Constraints\" must be a list with expressions of form ",
"0<term or 0=term; it may be empty")]; RETURN(ListOfErrors) end if; if not 
type(Sys["Variables"],list(name = EvalsToFloat .. EvalsToFloat)) then 
ListOfErrors := [op(ListOfErrors), "entry \"Variables\" must be a list with ex\
pressions of type name= EvalsToFloat..EvalsToFloat"]; RETURN(ListOfErrors) end
if; if not type(Sys["Parameters"],list(name = EvalsToFloat)) then ListOfErrors
:= [op(ListOfErrors), 
"entry \"Parameters\" must be a list of expressions of form name= EvalsToFloat\
"] end if; if not type(Sys["CostFunction"],{[], [term]}) then ListOfErrors :=
[op(ListOfErrors), 
"entry \"CostFunction\" must be a list with one entry of the form [term]"] end
if; if not ListOfErrors = [] then return ListOfErrors end if; if 1 < nargs 
then if args[2] = ('long') or args[2] = ('all') or args[2] = ('strict') then 
ListOfErrors := listOfErrorsInNlpPart2(Sys) end if end if; return ListOfErrors
end proc; noNameConflictsInUnitsInDAESys := proc (aDAESys) local VarNames, 
UnitNames, Conflicts; if not member("Units",map(op,[indices(aDAESys)])) then 
error "entry Units must exist in DAESys" end if; UnitNames := map(ModelPack:-
GetListOfIndetsIn,aDAESys["Units"]); UnitNames := map(op,UnitNames); UnitNames
:= convert(UnitNames,set); VarNames := [op(map(lhs,aDAESys["Parameters"])), op
(aDAESys["DynVars"]), op(aDAESys["AlgVars"]), op(map(lhs,aDAESys["ExplicitAEs"
]))]; VarNames := convert(VarNames,set); Conflicts := `intersect`(UnitNames,
VarNames); if not Conflicts = {} then error 
"conflicts in unit names and variables names for %1", Conflicts end if; return
true end proc; noNameConflictsInUnitsInDTASys := proc (aDTASys) local VarNames
, UnitNames, Conflicts; if not member("Units",map(op,[indices(aDTASys)])) then
error "entry Units must exist in DTASys" end if; UnitNames := map(Aux:-
ListOperations:-getListOfIndetsIn,aDTASys["Units"]); UnitNames := map(op,
UnitNames); UnitNames := convert(UnitNames,set); VarNames := [op(map(lhs,
aDTASys["Parameters"])), op(aDTASys["DynVars"]), op(aDTASys["AlgVars"]), op(
map(lhs,aDTASys["ExplicitAlgEqns"]))]; VarNames := convert(VarNames,set); 
Conflicts := `intersect`(UnitNames,VarNames); if not Conflicts = {} then error
"conflicts in unit names and variables names for %1", Conflicts end if; return
true end proc; subsExplicitAEsIntoDAESys := proc (aDAESys::DAESys) local i1, 
item, ListOfUnaffectedEntries, FormerODEs, FormerExplicitAEs, FormerAEs, 
PositionODEs, PositionExplicitAEs, NewODEs, NewExplicitAEs, NewDAESys, 
LHSNewDAESys; NewDAESys := table(); ListOfUnaffectedEntries := map(op,[indices
(aDAESys)]); ListOfUnaffectedEntries := Aux:-ListOperations:-
removeItemFromList(["ExplicitAEs", "ODEs", "AEs"],ListOfUnaffectedEntries); 
for item in ListOfUnaffectedEntries while true do NewDAESys[item] := aDAESys[
item] end do; NewDAESys["ExplicitAEs"] := []; NewDAESys["ODEs"] := Aux:-
ListOperations:-subsEqnListIntoEqn(aDAESys["ExplicitAEs"],aDAESys["ODEs"]); 
NewDAESys["AEs"] := Aux:-ListOperations:-subsEqnListIntoEqn(aDAESys[
"ExplicitAEs"],aDAESys["AEs"]); RETURN(eval(NewDAESys)) end proc; 
subsExplicitAlgEqnsIntoDTASys := proc (aDTASys::table) local i1, item, 
ListOfUnaffectedEntries, FormerDynEqns, FormerExplicitAlgEqns, FormerAlgEqns,
PositionDynEqns, PositionExplicitAlgEqns, NewDynEqns, NewExplicitAlgEqns, 
NewDTASys, LHSNewDTASys; NewDTASys := table(); ListOfUnaffectedEntries := map(
op,[indices(aDTASys)]); ListOfUnaffectedEntries := Aux:-ListOperations:-
removeItemFromList(["ExplicitAlgEqns", "DynEqns", "AlgEqns"],
ListOfUnaffectedEntries); for item in ListOfUnaffectedEntries while true do 
NewDTASys[item] := aDTASys[item] end do; NewDTASys["ExplicitAlgEqns"] := []; 
NewDTASys["DynEqns"] := Aux:-ListOperations:-subsEqnListIntoEqn(aDTASys[
"ExplicitAlgEqns"],aDTASys["DynEqns"]); NewDTASys["AlgEqns"] := Aux:-
ListOperations:-subsEqnListIntoEqn(aDTASys["ExplicitAlgEqns"],aDTASys[
"AlgEqns"]); RETURN(eval(NewDTASys)) end proc; end module; TensProd := module
() export Tijk_xj, Tijk_Akl, Trans_Tijk_Tikj, Sij_Tjk, Tij_xj, Aij_xj, xi_Tij,
xi_Aij, xi_Tij_yj, xi_Aij_yj; Tijk_xj := proc (Amat::('array')(3), xVec::{
Vector, list}) local x, ListOfArrayRanges, LengthOfVec, Result, i1, i2, i3; if
not type(xVec,Vector) then x := convert(xVec,Vector) else x := xVec end if; 
ListOfArrayRanges := op(2,eval(Amat)); LengthOfVec := nops(convert(x,list)); 
if not op(2,ListOfArrayRanges[2]) = LengthOfVec then error 
"incompatible dimensions of 1st and 2nd argument" end if; Result := array(
ListOfArrayRanges[1],ListOfArrayRanges[3]); for i1 to op(2,ListOfArrayRanges[1
]) do for i3 to op(2,ListOfArrayRanges[3]) do Result[i1,i3] := add(Amat[i1,i2,
i3]*x[i2],i2 = 1 .. LengthOfVec) end do end do; return eval(Result) end proc;
Tijk_Akl := proc (Ten3::('array')(3), Amat::{Matrix, ('array')(2)}) local 
ListOfTen3Ranges, ListOfAmatRanges, Result, i1, i2, i3, i4; ListOfTen3Ranges 
:= op(2,eval(Ten3)); ListOfAmatRanges := op(2,eval(Amat)); if not op(2,
ListOfAmatRanges[1]) = op(2,ListOfTen3Ranges[3]) then error 
"incompatible dimensions of 1st and 2nd argument" end if; Result := array(
ListOfTen3Ranges[1],ListOfTen3Ranges[2],ListOfAmatRanges[2]); for i1 to op(2,
ListOfTen3Ranges[1]) do for i2 to op(2,ListOfTen3Ranges[2]) do for i4 to op(2,
ListOfAmatRanges[2]) do Result[i1,i2,i4] := add(Ten3[i1,i2,i3]*Amat[i3,i4],i3
= 1 .. op(2,ListOfAmatRanges[1])) end do end do end do; return eval(Result) 
end proc; Trans_Tijk_Tikj := proc (Ten3::('array')(3)) local ListOfTen3Ranges,
Result, i1, i2, i3; ListOfTen3Ranges := op(2,eval(Ten3)); Result := array(
ListOfTen3Ranges[1],ListOfTen3Ranges[3],ListOfTen3Ranges[2]); for i1 to op(2,
ListOfTen3Ranges[1]) do for i2 to op(2,ListOfTen3Ranges[2]) do for i3 to op(2,
ListOfTen3Ranges[3]) do Result[i1,i3,i2] := Ten3[i1,i2,i3] end do end do end 
do; return eval(Result) end proc; Sij_Tjk := proc (Amat::{Matrix, ('array')(2)
}, Bmat::{Matrix, ('array')(1), ('array')(2)}) local i, j, k, A, B, NumRowsA,
NumColsA, NumRowsB, NumColsB, Sij_tjkInstance; if not type(Amat,Matrix) then A
:= convert(Amat,Matrix) else A := Amat end if; if not type(Bmat,Matrix) then B
:= convert(Bmat,Matrix) else B := Bmat end if; NumRowsA := nops(convert(A,
listlist)); NumColsA := nops(convert(A,listlist)[1]); NumRowsB := nops(convert
(B,listlist)); NumColsB := nops(convert(B,listlist)[1]); if not NumColsA = 
NumRowsB then error "incompatible dimensions of 1st and 2nd argument" end if;
Sij_tjkInstance := array(1 .. NumRowsA,1 .. NumColsB); for i to NumRowsA do 
for j to NumColsB do Sij_tjkInstance[i,j] := add(A[i,k]*B[k,j],k = 1 .. 
NumColsA) end do end do; return eval(Sij_tjkInstance) end proc; Tij_xj := proc
(Amat::{Matrix, array}, xVec::{Vector, list}) local i1, x, A, NumColsA, 
NumRowsA, ResultInArray, ResultInList; if not type(xVec,Vector) then x := 
convert(xVec,Vector) else x := xVec end if; if not type(Amat,Matrix) then A :=
convert(Amat,Matrix) else A := Amat end if; NumColsA := nops(convert(A,
listlist)[1]); if not NumColsA = nops(convert(x,list)) then error 
"incompatible dimensions of 1st and 2nd argument" end if; ResultInArray := 
Sij_Tjk(A,x); NumRowsA := nops(convert(A,listlist)); ResultInList := [seq(
ResultInArray[i1,1],i1 = 1 .. NumRowsA)]; return ResultInList end proc; Aij_xj
:= proc (Amat::{Matrix, array}, xVec::{Vector, list}) local A, NumRows, 
NumCols, result; if not type(Amat,Matrix) then A := convert(Amat,Matrix) else
A := Amat end if; NumRows := nops(convert(A,listlist)); NumCols := nops(
convert(A,listlist)[1]); if not NumRows = NumCols then error 
"1st  argument is not square Matrix." end if; result := Tij_xj(A,xVec); return
eval(result) end proc; xi_Tij := proc (xVec::{Vector, list}, Amat::{Matrix, ('
array')(2)}) local j1, x, A, NumRowsA, NumColsA, y, ResultInArray, 
ResultInList; if not type(xVec,list) then x := convert(xVec,list) else x := 
xVec end if; if not type(Amat,Matrix) then A := convert(Amat,Matrix) else A :=
Amat end if; NumRowsA := nops(convert(A,listlist)); NumColsA := nops(convert(A
,listlist)[1]); if not NumRowsA = nops(x) then error 
"incompatible dimensions of 1st and 2nd argument" end if; y := array([x]); 
ResultInArray := Sij_Tjk(y,A); ResultInList := [seq(ResultInArray[1,j1],j1 = 1
.. NumColsA)]; return ResultInList end proc; xi_Aij := proc (xVec::{Vector, 
list}, Amat::{Matrix, array}) local A, NumRows, NumCols, result; if not type(
Amat,Matrix) then A := convert(Amat,Matrix) else A := Amat end if; NumRows :=
nops(convert(A,listlist)); NumCols := nops(convert(A,listlist)[1]); if not 
NumRows = NumCols then error "2nd  argument is not square Matrix." end if; 
result := xi_Tij(xVec,A); return eval(result) end proc; xi_Tij_yj := proc (
xVec::{Vector, list}, Amat::{Matrix, ('array')(2)}, yVec::{Vector, list}) 
local x, y, A, NumRowsA, NumColsA, media, result; if not type(xVec,list) then
x := convert(xVec,list) else x := xVec end if; if not type(Amat,Matrix) then A
:= convert(Amat,Matrix) else A := Amat end if; if not type(yVec,Vector) then y
:= convert(yVec,Vector) else y := yVec end if; NumRowsA := nops(convert(A,
listlist)); if not NumRowsA = nops(x) then error 
"incompatible dimensions of 1st and 2nd argument" end if; NumColsA := nops(
convert(A,listlist)[1]); if not NumColsA = nops(convert(y,list)) then error 
"incompatible dimensions of 2nd and 3rd argument" end if; media := array([
xi_Tij(x,A)]); result := Sij_Tjk(media,y)[1,1]; return eval(result) end proc;
xi_Aij_yj := proc (xVec::{Vector, list}, Amat::{Matrix, ('array')(2)}, yVec::{
Vector, list}) local A, NumRows, NumCols, result; if not type(Amat,Matrix) 
then A := convert(Amat,Matrix) else A := Amat end if; NumRows := nops(convert(
A,listlist)); NumCols := nops(convert(A,listlist)[1]); if not NumRows = 
NumCols then error "1st  argument is not square Matrix." end if; result := 
xi_Tij_yj(xVec,A,yVec); return eval(result) end proc; end module; 
TransferNlpToGams := module () export convertMapleExprToGamsExpr, 
createFileWithTextForGams, createInitialPointTextForGams, createTextForGams, 
renameIndexedNameForGams, renameIndexedParamsAndVarsInNLP; 
convertMapleExprToGamsExpr := proc (mapleExpr) local cCodeExprString, 
indexOfEqSign, isMatchRegExpr, suffNumForComma, subStrForSubstitution, 
suffNumForPower, SubStr1, SubStr2, SubStr3, numForCommaList, numForBracket, 
powerExpr, powerExprNumeric, idexesForPower, SubStr1New, srtUnderPower, i; 
cCodeExprString := CodeGeneration[C](mapleExpr,resultname = "tempCodeCGenName"
,output = string); indexOfEqSign := StringTools[FirstFromLeft]("=",
cCodeExprString); cCodeExprString := StringTools[Delete](cCodeExprString,1 ..
indexOfEqSign); cCodeExprString := StringTools[SubstituteAll](cCodeExprString,
";",""); cCodeExprString := StringTools[Trim](cCodeExprString); 
cCodeExprString := StringTools[SubstituteAll](cCodeExprString,"(double) ","");
cCodeExprString := StringTools[SubstituteAll](cCodeExprString,"(int) ",""); 
cCodeExprString := StringTools[SubstituteAll](cCodeExprString,"pow","power");
isMatchRegExpr := StringTools[RegMatch]("power( )*(.*, -.*)",cCodeExprString);
 while isMatchRegExpr do suffNumForComma := StringTools[Search](", -",
cCodeExprString); subStrForSubstitution := StringTools[SubString](
cCodeExprString,1 .. suffNumForComma-1); suffNumForPower := [StringTools[
SearchAll]("power",subStrForSubstitution)]; SubStr1 := StringTools[SubString](
subStrForSubstitution,1 .. suffNumForPower[-1]-1); SubStr2 := StringTools[
SubString](subStrForSubstitution,suffNumForPower[-1]+5 .. StringTools[Length](
subStrForSubstitution)); SubStr3 := StringTools[SubString](cCodeExprString,
suffNumForComma+3 .. StringTools[Length](cCodeExprString)); cCodeExprString :=
StringTools[Join]([SubStr1, "1/power", SubStr2, ", ", SubStr3],""); 
isMatchRegExpr := StringTools[RegMatch]("power( )*(.*, -.*)",cCodeExprString)
end do; isMatchRegExpr := StringTools[RegMatch]("power( )*(.*, 0.1e1 / 0.2e1)"
,cCodeExprString);  while isMatchRegExpr do suffNumForComma := StringTools[
Search](", 0.1e1 / 0.2e1",cCodeExprString); subStrForSubstitution := 
StringTools[SubString](cCodeExprString,1 .. suffNumForComma-1); 
suffNumForPower := [StringTools[SearchAll]("power",subStrForSubstitution)]; 
SubStr1 := StringTools[SubString](subStrForSubstitution,1 .. suffNumForPower[-\
1]-1); SubStr2 := StringTools[SubString](subStrForSubstitution,suffNumForPower
[-1]+5 .. StringTools[Length](subStrForSubstitution)); SubStr3 := StringTools[
SubString](cCodeExprString,suffNumForComma+15 .. StringTools[Length](
cCodeExprString)); cCodeExprString := StringTools[Join]([SubStr1, "sqrt", 
SubStr2, SubStr3],""); isMatchRegExpr := StringTools[RegMatch](
"power( )*(.*, 0.1e1 / 0.2e1)",cCodeExprString) end do; isMatchRegExpr := 
StringTools[RegMatch]("power( )*(.*, 0.2e1)",cCodeExprString);  while 
isMatchRegExpr do suffNumForComma := StringTools[Search](", 0.2e1",
cCodeExprString); subStrForSubstitution := StringTools[SubString](
cCodeExprString,1 .. suffNumForComma-1); suffNumForPower := [StringTools[
SearchAll]("power",subStrForSubstitution)]; SubStr1 := StringTools[SubString](
subStrForSubstitution,1 .. suffNumForPower[-1]-1); SubStr2 := StringTools[
SubString](subStrForSubstitution,suffNumForPower[-1]+5 .. StringTools[Length](
subStrForSubstitution)); SubStr3 := StringTools[SubString](cCodeExprString,
suffNumForComma+7 .. StringTools[Length](cCodeExprString)); cCodeExprString :=
StringTools[Join]([SubStr1, "sqr", SubStr2, SubStr3],""); isMatchRegExpr := 
StringTools[RegMatch]("power( )*(.*, 0.2e1)",cCodeExprString) end do; 
numForCommaList := [StringTools[SearchAll](", ",cCodeExprString)]; for i to 
nops(numForCommaList) do numForCommaList := [StringTools[SearchAll](",",
cCodeExprString)]; SubStr1 := StringTools[SubString](cCodeExprString,1 .. 
numForCommaList[i]); subStrForSubstitution := StringTools[SubString](
cCodeExprString,numForCommaList[i]+2 .. StringTools[Length](cCodeExprString));
numForBracket := StringTools[Search](")",subStrForSubstitution); powerExpr :=
StringTools[SubString](subStrForSubstitution,1 .. numForBracket-1); 
powerExprNumeric := parse(powerExpr); SubStr3 := StringTools[SubString](
subStrForSubstitution,numForBracket .. StringTools[Length](
subStrForSubstitution)); if trunc(powerExprNumeric)-powerExprNumeric = 0 then
SubStr2 := convert(trunc(powerExprNumeric),string); cCodeExprString := 
StringTools[Join]([SubStr1, SubStr2, SubStr3],"") else SubStr2 := convert(
evalf(powerExprNumeric,5),string); idexesForPower := [StringTools[SearchAll](
"power",SubStr1)]; SubStr1New := StringTools[SubString](SubStr1,1 .. 
idexesForPower[-1]-1); srtUnderPower := StringTools[SubString](SubStr1,
idexesForPower[-1]+6 .. StringTools[Length](SubStr1)-1); cCodeExprString := 
StringTools[Join]([SubStr1New, "exp(", "(", SubStr2, ")", "*(", "log(", 
srtUnderPower, ")", ")", SubStr3],"") end if end do; return cCodeExprString 
end proc; createFileWithTextForGams := proc (aNLP::NLP, filePathGAMS::string)
local fd, textForGAMS, fileName, modelName; if 2 < nargs then modelName := 
convert(args[3],string) else fileName := FileTools[Filename](filePathGAMS); 
modelName := StringTools[Split](fileName,".")[1] end if; textForGAMS := Aux:-
TransferNlpToGams:-createTextForGams(aNLP,modelName); fd := fopen(filePathGAMS
,WRITE,TEXT); fprintf(fd,textForGAMS); fclose(fd); return "file was created" 
end proc; createInitialPointTextForGams := proc (initPoint::{list(name = 
EvalsToFloat)}) local counter1, varForGams, initPointGamsText; 
initPointGamsText := ""; for counter1 to nops(initPoint) do varForGams := 
renameIndexedNameForGams(lhs(initPoint[counter1])); initPointGamsText := 
StringTools[Join]([initPointGamsText, convert(varForGams,string), ".l=", 
convert(rhs(initPoint[counter1]),string), ";\n"],"") end do; return 
initPointGamsText end proc; createTextForGams := proc (aNLPpar::NLP, modelName
::string) local finalText, mapleParameters, i, textGamsParameters, 
newLineTextGamsParameters, mapleVariables, textGamsVariables, 
newLineTextGamsVariables, newLineTextLowerBound, newLineTextUpperBound, 
mapleCurrentVar, mapleCurrentLowerBound, mapleCurrentUpperBound, 
mapleIntegerVariables, textGamsIntegerVariables, mapleCurrentIntVar, 
mapleCurrentIntLowerBound, mapleCurrentIntUpperBound, mapleBinaryVariables, 
textGamsBinaryVariables, newLineTextGamsBinaryVariables, 
newLineTextGamsIntVariables, newLineTextLowerIntBound, 
newLineTextUpperIntBound, newLineTextLowerBinaryBound, 
newLineTextUpperBinaryBound, mapleConstraints, textGamsDeclarationEquations, 
currentGamsConstraintLine, textGamsConstraints, signOfConstraint, mapleObjFunc
, textGamsForObjFuncConstraint, iflessConstrRhs, iflessOrEqConstrRhs, 
ifGraterConstrRhs, ifGraterOrEqConstrRhs, iflessConstrLhs, textGamsModelName,
textGamsSolverOption1, textGamsSolverOption2, textGamsSolveCommand, 
textGamsDisplayCommand, ifHasIntVariables, solverType, aNLP; aNLP := Aux:-
TransferNlpToGams:-renameIndexedParamsAndVarsInNLP(aNLPpar); mapleParameters 
:= aNLP["Parameters"]; textGamsParameters := ""; for i to nops(mapleParameters
) do newLineTextGamsParameters := StringTools[Join](["Scalar ", convert(lhs(
mapleParameters[i]),string), " /", convert(rhs(mapleParameters[i]),string), 
"/;"],""); textGamsParameters := StringTools[Join]([textGamsParameters, 
newLineTextGamsParameters, "\n"],"") end do; mapleVariables := aNLP[
"Variables"]; textGamsVariables := ""; for i to nops(mapleVariables) do 
mapleCurrentVar := lhs(mapleVariables[i]); mapleCurrentLowerBound := lhs(rhs(
mapleVariables[i])); mapleCurrentUpperBound := rhs(rhs(mapleVariables[i])); 
newLineTextGamsVariables := StringTools[Join](["Variable ", convert(
mapleCurrentVar,string), ";\n"],""); newLineTextLowerBound := ""; 
newLineTextUpperBound := ""; if -infinity < mapleCurrentLowerBound then 
newLineTextLowerBound := StringTools[Join](["   ", convert(mapleCurrentVar,
string), ".lo=", convert(mapleCurrentLowerBound,string), ";\n"],"") end if; if
mapleCurrentUpperBound < infinity then newLineTextUpperBound := StringTools[
Join](["   ", convert(mapleCurrentVar,string), ".up=", convert(
mapleCurrentUpperBound,string), ";\n"],"") end if; textGamsVariables := 
StringTools[Join]([textGamsVariables, newLineTextGamsVariables, 
newLineTextLowerBound, newLineTextUpperBound, "\n"],"") end do; 
textGamsVariables := StringTools[Join]([textGamsVariables, 
"Variable SlackVariableForMinimizing;\n", 
"SlackVariableForMinimizing.lo=-100000;\n", 
"SlackVariableForMinimizing.up=100000;\n\n"],""); mapleIntegerVariables := 
aNLP["IntegerVariables"]; textGamsIntegerVariables := ""; ifHasIntVariables :=
false; if type(mapleIntegerVariables,list(name = EvalsToFloat .. EvalsToFloat)
) then ifHasIntVariables := true; for i to nops(mapleIntegerVariables) do 
mapleCurrentIntVar := lhs(mapleIntegerVariables[i]); mapleCurrentIntLowerBound
:= lhs(rhs(mapleIntegerVariables[i])); mapleCurrentIntUpperBound := rhs(rhs(
mapleIntegerVariables[i])); newLineTextGamsIntVariables := StringTools[Join]([
"Integer variable ", convert(mapleCurrentIntVar,string), ";\n"],""); 
newLineTextLowerIntBound := ""; newLineTextUpperIntBound := ""; if 0 < 
mapleCurrentIntLowerBound then newLineTextLowerIntBound := StringTools[Join]([
"   ", convert(mapleCurrentIntVar,string), ".lo=", convert(
mapleCurrentIntLowerBound,string), ";\n"],"") end if; if 
mapleCurrentIntUpperBound < infinity then newLineTextUpperIntBound := 
StringTools[Join](["   ", convert(mapleCurrentIntVar,string), ".up=", convert(
mapleCurrentIntUpperBound,string), ";\n"],"") end if; textGamsIntegerVariables
:= StringTools[Join]([textGamsIntegerVariables, newLineTextGamsIntVariables, 
newLineTextLowerIntBound, newLineTextUpperIntBound, "\n"],"") end do end if; 
mapleBinaryVariables := aNLP["BinaryVariables"]; textGamsBinaryVariables := ""
; if type(mapleBinaryVariables,list(name)) then ifHasIntVariables := true; for
i to nops(mapleBinaryVariables) do newLineTextGamsBinaryVariables := 
StringTools[Join](["Binary variable ", convert(mapleBinaryVariables[i],string)
, ";\n"],""); newLineTextLowerBinaryBound := StringTools[Join](["   ", convert
(mapleBinaryVariables[i],string), ".lo=0", ";\n"],""); 
newLineTextUpperBinaryBound := StringTools[Join](["   ", convert(
mapleBinaryVariables[i],string), ".up=1", ";\n"],""); textGamsBinaryVariables
:= StringTools[Join]([textGamsBinaryVariables, newLineTextGamsBinaryVariables,
newLineTextLowerBinaryBound, newLineTextUpperBinaryBound, "\n"],"") end do end
if; mapleConstraints := aNLP["Constraints"]; textGamsDeclarationEquations := 
"Equations \n"; for i to nops(mapleConstraints) do 
textGamsDeclarationEquations := StringTools[Join]([
textGamsDeclarationEquations, "   EquationN", convert(i,string), ",\n"],"") 
end do; textGamsDeclarationEquations := StringTools[Join]([
textGamsDeclarationEquations, "   ObjectiveFunctionEquation;\n"],""); 
textGamsConstraints := ""; for i to nops(mapleConstraints) do signOfConstraint
:= " =e= "; iflessConstrLhs := convert(mapleConstraints[i],string); 
iflessConstrRhs := convert(lhs(mapleConstraints[i]) < rhs(mapleConstraints[1])
,string); iflessOrEqConstrRhs := convert(lhs(mapleConstraints[i]) <= rhs(
mapleConstraints[1]),string); if StringTools[Compare](iflessConstrLhs,
iflessConstrRhs) or StringTools[Compare](iflessConstrLhs,iflessOrEqConstrRhs)
then signOfConstraint := " =l= " end if; ifGraterConstrRhs := convert(rhs(
mapleConstraints[1]) < lhs(mapleConstraints[i]),string); ifGraterOrEqConstrRhs
:= convert(rhs(mapleConstraints[1]) <= lhs(mapleConstraints[i]),string); if 
StringTools[Compare](iflessConstrLhs,ifGraterConstrRhs) or StringTools[Compare
](iflessConstrLhs,ifGraterOrEqConstrRhs) then signOfConstraint := " =g= " end
if; currentGamsConstraintLine := StringTools[Join](["EquationN", convert(i,
string), " .. ", Aux:-TransferNlpToGams:-convertMapleExprToGamsExpr(lhs(
mapleConstraints[i])), signOfConstraint, Aux:-TransferNlpToGams:-
convertMapleExprToGamsExpr(rhs(mapleConstraints[i])), ";"],""); 
textGamsConstraints := StringTools[Join]([textGamsConstraints, 
currentGamsConstraintLine, "\n"],"") end do; mapleObjFunc := aNLP[
"CostFunction"]; textGamsForObjFuncConstraint := StringTools[Join]([
"ObjectiveFunctionEquation", " .. ", "SlackVariableForMinimizing", " =e= ", 
Aux:-TransferNlpToGams:-convertMapleExprToGamsExpr(mapleObjFunc), ";\n"],"");
textGamsModelName := StringTools[Join](["Model ", modelName, " /all/;\n"],"");
textGamsSolverOption1 := "option nlp=SNOPT;\n"; textGamsSolverOption2 := 
"option minlp=SBB;\n"; solverType := "nlp"; if ifHasIntVariables = true then 
solverType := "minlp" end if; textGamsSolveCommand := StringTools[Join]([
"solve ", modelName, " using ", solverType, 
" minimizing SlackVariableForMinimizing;\n"],""); textGamsDisplayCommand := 
"display SlackVariableForMinimizing.l;"; finalText := StringTools[Join]([
textGamsParameters, "\n", textGamsVariables, textGamsIntegerVariables, 
textGamsBinaryVariables, "\n", textGamsDeclarationEquations, "\n", 
textGamsConstraints, textGamsForObjFuncConstraint, "\n", textGamsModelName, 
textGamsSolverOption1, textGamsSolverOption2, textGamsSolveCommand, 
textGamsDisplayCommand],""); RETURN(finalText) end proc; 
renameIndexedNameForGams := proc (varToRename::name) local varToRenameAsString
, newVarName; varToRenameAsString := convert(varToRename,string); 
varToRenameAsString := StringTools[SubstituteAll](varToRenameAsString,"[","_")
; varToRenameAsString := StringTools[SubstituteAll](varToRenameAsString,"]",""
); varToRenameAsString := StringTools[SubstituteAll](varToRenameAsString,",",
""); newVarName := convert(varToRenameAsString,symbol); RETURN(newVarName) end
proc; renameIndexedParamsAndVarsInNLP := proc (aNLP::NLP) local 
mapleParameters, mapleVariables, mapleIntegerVariables, mapleBinaryVariables,
hasIntegerVariables, hasBinaryVariables, varsAndPars, i, newVarName, 
listOfNaamesForSubstitution, newParam, varToRename, newVars, newIntVars, 
newBinVars, newContrs, newCostFun, newAEs, newLinContrs; mapleParameters := 
aNLP["Parameters"]; mapleVariables := aNLP["Variables"]; mapleIntegerVariables
:= aNLP["IntegerVariables"]; hasIntegerVariables := false; 
mapleBinaryVariables := aNLP["BinaryVariables"]; hasBinaryVariables := false;
varsAndPars := [op(map(lhs,mapleParameters)), op(map(lhs,mapleVariables))]; if
type(mapleIntegerVariables,list(name = EvalsToFloat .. EvalsToFloat)) then 
hasIntegerVariables := true; varsAndPars := [op(varsAndPars), op(map(lhs,
mapleIntegerVariables))] end if; if type(mapleBinaryVariables,list(name)) then
hasBinaryVariables := true; varsAndPars := [op(varsAndPars), op(
mapleBinaryVariables)] end if; listOfNaamesForSubstitution := []; for i to 
nops(varsAndPars) do if type(varsAndPars[i],indexed) then varToRename := 
varsAndPars[i]; newVarName := renameIndexedNameForGams(varToRename); 
listOfNaamesForSubstitution := [op(listOfNaamesForSubstitution), varToRename =
newVarName] end if end do; newParam := subs(listOfNaamesForSubstitution,
mapleParameters); newVars := subs(listOfNaamesForSubstitution,mapleVariables);
if hasIntegerVariables = true then newIntVars := subs(
listOfNaamesForSubstitution,mapleIntegerVariables) else newIntVars := [] end 
if; if hasBinaryVariables = true then newBinVars := subs(
listOfNaamesForSubstitution,mapleBinaryVariables) else newBinVars := [] end if
; newContrs := subs(listOfNaamesForSubstitution,aNLP["Constraints"]); 
newCostFun := subs(listOfNaamesForSubstitution,aNLP["CostFunction"]); if type(
anNLP["ExplicitAEs"],list) then newAEs := subs(listOfNaamesForSubstitution,
aNLP["ExplicitAEs"]) else newAEs := [] end if; if type(anNLP[
"LinearConstraints"],list) then newLinContrs := subs(
listOfNaamesForSubstitution,aNLP["LinearConstraints"]) else newLinContrs := []
end if; RETURN(table(["Parameters" = newParam, "Variables" = newVars, 
"IntegerVariables" = newIntVars, "BinaryVariables" = newBinVars, "Constraints"
= newContrs, "CostFunction" = newCostFun, "ExplicitAEs" = newAEs, 
"LinearConstraints" = newLinContrs])) end proc; end module; end module;

  Aux:-init();


