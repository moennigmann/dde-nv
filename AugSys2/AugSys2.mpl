AugSys2 := module () export ConstraintInitializer, Discretization, 
FlipAndFoldNV, HopfNV, NeimarkSackerNV, SaddleNodeNV; ConstraintInitializer :=
module () local augSysInModule, isAugSysCreatedInModule, 
listOfNVParamsInModule, listOfDynVarsInModule, isNSType, isFlipOrFoldType, 
isSaddleNodeType, isHopfType; export calcClosestCriticalPoint, calcOptPoint, 
createAugSys, CriticalPoint, FixPoint, getAugSys, getScaledDTASys, NLP, 
runOptForListStartPoints; calcClosestCriticalPoint := proc () local anNLP, 
ListOfErrorsInAugSysNLP, NPSOLproc, StrartingCriticalPoint, Inform, Objf, 
OptPoint, isCalcClosestCPSucceeded; isCalcClosestCPSucceeded := false; AugSys2
:-ConstraintInitializer:-CriticalPoint:-setCalcClosestCPSucceeded(
isCalcClosestCPSucceeded); anNLP := AugSys2:-ConstraintInitializer:-NLP:-
getNLP(); ListOfErrorsInAugSysNLP := Aux:-SystemClasses:-listOfErrorsInNLP(
anNLP,'strict'); if ListOfErrorsInAugSysNLP <> true then error 
"NLP is not defined correct", ListOfErrorsInAugSysNLP end if; if Aux:-
FileOperations:-dirExists("tmp") then error 
"temporary directory ./tmp needed to run test already exists" else mkdir("tmp"
) end if; StrartingCriticalPoint := AugSys2:-ConstraintInitializer:-
CriticalPoint:-getStartingCriticalPoint(); if not type(StrartingCriticalPoint,
list(EvalsToFloat)) then error "Starting critical point is not defined" end if
; printf("\n"); printf("running Optimization for NLP\n"); NPSOLproc := NPSOL:-
CreateInstance(anNLP,"./tmp"); NPSOLproc:-setXVEC(StrartingCriticalPoint); 
NPSOLproc:-runOpt(); Inform := NPSOLproc:-getInform(); if Inform = 0 then 
isCalcClosestCPSucceeded := true; AugSys2:-ConstraintInitializer:-
CriticalPoint:-setCalcClosestCPSucceeded(isCalcClosestCPSucceeded); Objf := 
NPSOLproc:-getObjf(); OptPoint := NPSOLproc:-getXVEC(); AugSys2:-
ConstraintInitializer:-CriticalPoint:-setCalcClosestCPObjf(Objf); AugSys2:-
ConstraintInitializer:-CriticalPoint:-setCalcClosestCPOptPoint(convert(
OptPoint,list)) end if; return isCalcClosestCPSucceeded end proc; calcOptPoint
:= proc () local anNLP, ListOfErrorsInAugSysNLP, NPSOLproc, 
StrartingCriticalPoint, Inform, Objf, OptPoint, isCalcOPSucceeded; 
isCalcOPSucceeded := false; AugSys2:-ConstraintInitializer:-CriticalPoint:-
setCalcOPSucceeded(isCalcOPSucceeded); anNLP := AugSys2:-ConstraintInitializer
:-NLP:-getNLP(); ListOfErrorsInAugSysNLP := Aux:-SystemClasses:-
listOfErrorsInNLP(anNLP,'strict'); if ListOfErrorsInAugSysNLP <> true then 
error "NLP is not defined correct", ListOfErrorsInAugSysNLP end if; if Aux:-
FileOperations:-dirExists("tmp2") then error 
"temporary directory ./tmp2 needed to run test already exists" else mkdir(
"tmp2") end if; StrartingCriticalPoint := AugSys2:-ConstraintInitializer:-
CriticalPoint:-getStartingCriticalPoint(); if not type(StrartingCriticalPoint,
list(EvalsToFloat)) then error "Starting critical point is not defined" end if
; printf("\n"); printf("running Optimization for NLP\n"); NPSOLproc := NPSOL:-
CreateInstance(anNLP,"./tmp2"); NPSOLproc:-setXVEC(StrartingCriticalPoint); 
NPSOLproc:-runOpt(); Inform := NPSOLproc:-getInform(); if Inform = 0 then 
isCalcOPSucceeded := true; AugSys2:-ConstraintInitializer:-CriticalPoint:-
setCalcOPSucceeded(isCalcOPSucceeded); Objf := NPSOLproc:-getObjf(); OptPoint
:= NPSOLproc:-getXVEC(); AugSys2:-ConstraintInitializer:-CriticalPoint:-
setCalcOPObjf(Objf); AugSys2:-ConstraintInitializer:-CriticalPoint:-
setCalcOPOptPoint(convert(OptPoint,list)) end if; return isCalcOPSucceeded end
proc; createAugSys := proc (aSys::{DAESys, DTASys}, listOfNVparams::list(name)
, keyword::string) local augSys, radiusOfEig; if keyword = "Neimark-Sacker" or
keyword = "NeimarkSacker" or keyword = "NS" or keyword = "Flip" or keyword = 
"flip" or keyword = "Fold" or keyword = "fold" then if 3 < nargs then 
radiusOfEig := args[4]; if not (0 <= radiusOfEig and radiusOfEig <= 1) then 
error "Radius of circle where eigenvalues have to lie (forth input) have to be\
 between 0 to 1" end if else radiusOfEig := 1 end if else if keyword = 
"SaddleNode" or keyword = "saddleNode" or keyword = "Saddle" or keyword = 
"saddle" or keyword = "Hopf" or keyword = "hopf" then if 3 < nargs then 
radiusOfEig := args[4]; if 0 < radiusOfEig then error 
"Bound on a real part of eigenvalue should be negative" end if else 
radiusOfEig := 0 end if else error "Name of bifurcation (third parameter) is w\
rong defined. The keyword has to be one of the folowing strings \"NS\",\"flip\
\",\"fold\",\"saddle\" or \"hopf." end if end if; isNSType := false; 
isFlipOrFoldType := false; isSaddleNodeType = false; isHopfType = false; if 
keyword = "Neimark-Sacker" or keyword = "NeimarkSacker" or keyword = "NS" then
augSys := AugSys2:-NeimarkSackerNV:-CreateNSWithAlgEqnsNVSys(aSys,
listOfNVparams,radiusOfEig):-getSys(); isNSType := true else if keyword = 
"Flip" or keyword = "flip" then augSys := AugSys2:-FlipAndFoldNV:-
CreateFlipWithAlgEqnsNVSys(aSys,listOfNVparams,radiusOfEig):-getSys(); 
isFlipOrFoldType := true else if keyword = "Fold" or keyword = "fold" then 
augSys := AugSys2:-FlipAndFoldNV:-CreateFoldWithAlgEqnsNVSys(aSys,
listOfNVparams,radiusOfEig):-getSys(); isFlipOrFoldType := true else if 
keyword = "SaddleNode" or keyword = "saddleNode" or keyword = "Saddle" or 
keyword = "saddle" then augSys := AugSys2:-SaddleNodeNV:-CreateSaddleNodeNVSys
(aSys,listOfNVparams,radiusOfEig):-getSys(); isSaddleNodeType := true else if
keyword = "Hopf" or keyword = "hopf" then augSys := AugSys2:-HopfNV:-
CreateHopfNVSys(aSys,listOfNVparams,radiusOfEig):-getSys(); isHopfType := true
else error "The keyword has to be one of the folowing strings \"NS\",\"flip\",\
\"fold\",\"saddle\" or \"hopf\"" end if end if end if end if end if; 
augSysInModule := augSys; isAugSysCreatedInModule := true; 
listOfNVParamsInModule := listOfNVparams; listOfDynVarsInModule := aSys[
"DynVars"]; return eval(augSys) end proc; CriticalPoint := module () local 
criticalPointParametersInModule, criticalPointVariablesInModule, 
startingCriticalPointInModule, isCalcClosestCPSucceededInModule, 
calcClosestCPObjfInModule, calcClosestCPOptPointInModule, 
isCalcOPSucceededInModule, calcOPObjfInModule, calcOPOptPointInModule; export
getCalcClosestCPObjf, getCalcClosestCPOptPoint, getCalcOPObjf, 
getCalcOPOptPoint, getConstraints, getNVConstraints, getParameters, 
getStartingCriticalPoint, getVariables, isCalcClosestCPSucceeded, 
isCalcOPSucceeded, setCalcClosestCPObjf, setCalcClosestCPOptPoint, 
setCalcClosestCPSucceeded, setCalcOPObjf, setCalcOPOptPoint, 
setCalcOPSucceeded, setParameters, setStartingCPWithDistance, 
setStartingCriticalPoint, setVariables; getCalcClosestCPObjf := proc (V) 
return eval(calcClosestCPObjfInModule) end proc; getCalcClosestCPOptPoint := 
proc () return eval(calcClosestCPOptPointInModule) end proc; getCalcOPObjf :=
proc (V) return eval(calcOPObjfInModule) end proc; getCalcOPOptPoint := proc (
) return eval(calcOPOptPointInModule) end proc; getConstraints := proc () 
local numOfFixSystemEqns, augSysEqns, newEquations; if isAugSysCreatedInModule
<> true then error 
"first define augmented system with the procedure createAugSys" end if; 
numOfFixSystemEqns := nops(listOfDynVarsInModule); augSysEqns := 
augSysInModule["Equations"]; if isNSType = true then newEquations := [seq(
augSysEqns[i1],i1 = numOfFixSystemEqns+1 .. 3*numOfFixSystemEqns+2)] end if; 
if isFlipOrFoldType = true then newEquations := [seq(augSysEqns[i1],i1 = 
numOfFixSystemEqns+1 .. 2*numOfFixSystemEqns+1)] end if; if isSaddleNodeType =
true then newEquations := [seq(augSysEqns[i1],i1 = numOfFixSystemEqns+1 .. 2*
numOfFixSystemEqns+1)] end if; if isHopfType = true then newEquations := [seq(
augSysEqns[i1],i1 = numOfFixSystemEqns+1 .. 3*numOfFixSystemEqns+2)] end if; 
return eval(newEquations) end proc; getNVConstraints := proc () local 
numOfFixSystemEqns, augSysEqns, newEquations, numOfNVParams; if 
isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
numOfFixSystemEqns := nops(listOfDynVarsInModule); numOfNVParams := nops(
listOfNVParamsInModule); augSysEqns := augSysInModule["Equations"]; if 
isNSType = true then newEquations := [seq(augSysEqns[i1],i1 = 3*
numOfFixSystemEqns+3 .. 6*numOfFixSystemEqns+4+numOfNVParams)] end if; if 
isFlipOrFoldType = true then newEquations := [seq(augSysEqns[i1],i1 = 2*
numOfFixSystemEqns+2 .. 4*numOfFixSystemEqns+2+numOfNVParams)] end if; if 
isSaddleNodeType = true then newEquations := [seq(augSysEqns[i1],i1 = 2*
numOfFixSystemEqns+2 .. 2*numOfFixSystemEqns+1+numOfNVParams)] end if; if 
isHopfType = true then newEquations := [seq(augSysEqns[i1],i1 = 3*
numOfFixSystemEqns+3 .. 6*numOfFixSystemEqns+4+numOfNVParams)] end if; return
eval(newEquations) end proc; getParameters := proc () return eval(
criticalPointParametersInModule) end proc; getStartingCriticalPoint := proc ()
return eval(startingCriticalPointInModule) end proc; getVariables := proc () 
return eval(criticalPointVariablesInModule) end proc; isCalcClosestCPSucceeded
:= proc () return eval(isCalcClosestCPSucceededInModule) end proc; 
isCalcOPSucceeded := proc () return eval(isCalcOPSucceededInModule) end proc;
setCalcClosestCPObjf := proc (ValueObjf::EvalsToFloat) 
calcClosestCPObjfInModule := ValueObjf; return end proc; 
setCalcClosestCPOptPoint := proc (ValueOptPoint::{list(EvalsToFloat), list(
name = EvalsToFloat)}) local i, NewOptPoint, NamesOfVars, NLPVars; if type(
ValueOptPoint) = list(name = EvalsToFloat) then calcClosestCPOptPointInModule
:= ValueOptPoint else NLPVars := AugSys2:-ConstraintInitializer:-NLP:-
getVariables(); NamesOfVars := map(lhs,NLPVars); if nops(ValueOptPoint) <> 
nops(NamesOfVars) then error 
"Optimal Point has not the same order as \"Variables\" in NLP" end if; 
NewOptPoint := []; for i to nops(ValueOptPoint) do NewOptPoint := [op(
NewOptPoint), NamesOfVars[i] = ValueOptPoint[i]] end do; 
calcClosestCPOptPointInModule := eval(NewOptPoint) end if; return end proc; 
setCalcClosestCPSucceeded := proc (ValueCalcClosestCPSucceeded::truefalse) 
isCalcClosestCPSucceededInModule := ValueCalcClosestCPSucceeded; return end 
proc; setCalcOPObjf := proc (ValueObjf::EvalsToFloat) calcOPObjfInModule := 
ValueObjf; return end proc; setCalcOPOptPoint := proc (ValueOptPoint::{list(
EvalsToFloat), list(name = EvalsToFloat)}) local i, NewOptPoint, NamesOfVars,
NLPVars; if type(ValueOptPoint) = list(name = EvalsToFloat) then 
calcOPOptPointInModule := ValueOptPoint else NLPVars := AugSys2:-
ConstraintInitializer:-NLP:-getVariables(); NamesOfVars := map(lhs,NLPVars); 
if nops(ValueOptPoint) <> nops(NamesOfVars) then error 
"Optimal Point has not the same order as \"Variables\" in NLP" end if; 
NewOptPoint := []; for i to nops(ValueOptPoint) do NewOptPoint := [op(
NewOptPoint), NamesOfVars[i] = ValueOptPoint[i]] end do; 
calcOPOptPointInModule := eval(NewOptPoint) end if; return end proc; 
setCalcOPSucceeded := proc (ValueCalcOPSucceeded::truefalse) 
isCalcOPSucceededInModule := ValueCalcOPSucceeded; return end proc; 
setParameters := proc (NewPars::list(name = EvalsToFloat)) local NewParsNames,
ParNames, Pars, Missing, Obsolete; if isAugSysCreatedInModule <> true then 
error "first define augmented system with the procedure createAugSys" end if;
NewParsNames := map(lhs,NewPars); NewParsNames := [op(NewParsNames), op(
listOfNVParamsInModule)]; Pars := augSysInModule["Parameters"]; ParNames := 
map(lhs,Pars); Missing, Obsolete := Aux:-ListOperations:-
getMissingAndObsoleteNames(NewParsNames,ParNames); if not Missing = {} then 
error "assignments are missing for %1", Missing end if; 
criticalPointParametersInModule := NewPars; return end proc; 
setStartingCPWithDistance := proc (NewNamesOfVars::list(name), Distanse::
EvalsToFloat) local VariblesOfNLP, calcClosestPointSucceded, NumOfNVParams, 
ClosestCP, RVector, NormRVector, NVVars, i, ListofNewVars, NewVarValue, 
ListOfAllVars; if isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
VariblesOfNLP := AugSys2:-ConstraintInitializer:-NLP:-getVariables(); if not 
type(VariblesOfNLP,list(name = EvalsToFloat .. EvalsToFloat)) or nops(
VariblesOfNLP) = 0 then error "first define \"Variables\" in NLP" end if; 
calcClosestPointSucceded := AugSys2:-ConstraintInitializer:-CriticalPoint:-
isCalcClosestCPSucceeded(); if calcClosestPointSucceded <> true then error 
"Calculation of closest critical point on the boundary was not succeed" end if
; NumOfNVParams := nops(listOfNVParamsInModule); if nops(NewNamesOfVars) <> 
NumOfNVParams then error "Names of Optimized Variables have to be the same ord\
er as the list of NV Parameters in augmented system", nops(
listOfNVParamsInModule) end if; ClosestCP := AugSys2:-ConstraintInitializer:-
CriticalPoint:-getCalcClosestCPOptPoint(); RVector := Vector(NumOfNVParams,
symbol = r); RVector := subs(ClosestCP,RVector); NormRVector := LinearAlgebra[
Norm](RVector,2); NVVars := subs(ClosestCP,listOfNVParamsInModule); 
ListofNewVars := []; for i to NumOfNVParams do NewVarValue := NVVars[i]+
Distanse*RVector[i]/NormRVector; ListofNewVars := [op(ListofNewVars), 
NewNamesOfVars[i] = NewVarValue] end do; ListofNewVars := [op(ListofNewVars),
dist = abs(Distanse)]; ListOfAllVars := [op(ClosestCP), op(ListofNewVars)]; 
AugSys2:-ConstraintInitializer:-CriticalPoint:-setStartingCriticalPoint(
ListOfAllVars); return end proc; setStartingCriticalPoint := proc (NewVars::{
list(EvalsToFloat), list(name = EvalsToFloat)}) local VariblesOfNLP, Missing,
Obsolete, NamesOfVariblesOfNLP, NamesOfVariblesOfStartingPoint, i, j, 
OrderedCriticalPoint, NamesOfVariblesOfNLPString, i1, 
NamesOfVariblesOfStartingPointString; VariblesOfNLP := AugSys2:-
ConstraintInitializer:-NLP:-getVariables(); if not type(VariblesOfNLP,list(
name = EvalsToFloat .. EvalsToFloat)) or nops(VariblesOfNLP) = 0 then error 
"first define \"Variables\" in NLP" end if; if nops(NewVars) <> nops(
VariblesOfNLP) then error 
"Starting critical point must be the same order as number of variables in NLP"
, nops(VariblesOfNLP) end if; if type(NewVars,list(EvalsToFloat)) then 
startingCriticalPointInModule := NewVars end if; if type(NewVars,list(name = 
EvalsToFloat)) then NamesOfVariblesOfNLP := map(lhs,VariblesOfNLP); 
NamesOfVariblesOfStartingPoint := map(lhs,NewVars); NamesOfVariblesOfNLPString
:= []; NamesOfVariblesOfStartingPointString := []; for i1 to nops(
NamesOfVariblesOfNLP) do NamesOfVariblesOfNLPString := [op(
NamesOfVariblesOfNLPString), convert(NamesOfVariblesOfNLP[i1],string)]; 
NamesOfVariblesOfStartingPointString := [op(
NamesOfVariblesOfStartingPointString), convert(NamesOfVariblesOfStartingPoint[
i1],string)] end do; Missing, Obsolete := Aux:-ListOperations:-
getMissingAndObsoleteNames(NamesOfVariblesOfStartingPointString,
NamesOfVariblesOfNLPString); if not Missing = {} then error 
"assignments are missing for %1", Missing end if; OrderedCriticalPoint := 
array(1 .. nops(NamesOfVariblesOfNLP)); for i to nops(
NamesOfVariblesOfNLPString) do for j to nops(
NamesOfVariblesOfStartingPointString) do if NamesOfVariblesOfNLPString[i] = 
NamesOfVariblesOfStartingPointString[j] then OrderedCriticalPoint[i] := rhs(
NewVars[j]) end if end do end do; startingCriticalPointInModule := convert(
OrderedCriticalPoint,list) end if; return end proc; setVariables := proc (
NewVars::list(name = EvalsToFloat)) local NewVarsNames, Missing, Obsolete; if
isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
NewVarsNames := map(lhs,NewVars); Missing, Obsolete := Aux:-ListOperations:-
getMissingAndObsoleteNames(NewVarsNames,listOfDynVarsInModule); if not Missing
= {} then error "assignments are missing for %1", Missing end if; 
criticalPointVariablesInModule := NewVars; return end proc; end module; 
FixPoint := module () local fixPointParametersInModule, 
fixPointVariablesInModule, isFixPointParamsInModule; export getConstraints, 
getParameters, setParameters; getConstraints := proc () local 
numOfFixSystemEqns, augSysEqns, newEquations; if isAugSysCreatedInModule <> 
true then error 
"first define augmented system with the procedure createAugSys" end if; 
numOfFixSystemEqns := nops(listOfDynVarsInModule); augSysEqns := 
augSysInModule["Equations"]; newEquations := [seq(augSysEqns[i1],i1 = 1 .. 
numOfFixSystemEqns)]; return eval(newEquations) end proc; getParameters := 
proc () return eval(fixPointParametersInModule) end proc; setParameters := 
proc (NewPars::list(name = EvalsToFloat)) local NewParsNames, ParNames, Pars,
Missing, Obsolete; if isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
NewParsNames := map(lhs,NewPars); NewParsNames := [op(NewParsNames), op(
listOfNVParamsInModule)]; Pars := augSysInModule["Parameters"]; ParNames := 
map(lhs,Pars); Missing, Obsolete := Aux:-ListOperations:-
getMissingAndObsoleteNames(NewParsNames,ParNames); if not Missing = {} then 
error "assignments are missing for %1", Missing end if; 
fixPointParametersInModule := NewPars; isFixPointParamsInModule := true; 
return end proc; end module; getAugSys := proc () return eval(augSysInModule)
end proc; getScaledDTASys := proc (aSys::DTASys, listScaledParam::listlist) 
local numScaledParams, listSubsScaling, i1, aSysDynEqns, aSysExplicitAlgEqns,
aSysAlgEqns; numScaledParams := nops(listScaledParam); listSubsScaling := [];
for i1 to numScaledParams do listSubsScaling := [op(listSubsScaling), 
listScaledParam[i1][1] = listScaledParam[i1][1]*listScaledParam[i1][2]] end do
; aSysDynEqns := subs(listSubsScaling,aSys["DynEqns"]); aSysExplicitAlgEqns :=
subs(listSubsScaling,aSys["ExplicitAlgEqns"]); aSysAlgEqns := subs(
listSubsScaling,aSys["AlgEqns"]); return table(["ExplicitAlgEqns" = 
aSysExplicitAlgEqns, "DynEqns" = aSysDynEqns, "AlgEqns" = aSysAlgEqns, 
"DynVars" = aSys["DynVars"], "AlgVars" = aSys["AlgVars"], "Parameters" = aSys[
"Parameters"]]) end proc; NLP := module () local NLPinModule; export 
addConstraint, addConstraints, addCriticalPointConstraints, 
addDistanceConstraints, addFixPointConstraints, addGeneralConstraints, 
addLinearConstraint, addNVConstraintsForTighterApprBox, addNVConstraints, 
addVariable, addVariableWithoutRange, copyCriticalPointParamsToNLP, 
createNewNLP, getConstraints, getCostFunction, getExplicitAEs, 
getLinearConstraints, getNLP, getObjfunTemplatePath, getParameters, 
getVariables, setBinaryVariables, setConstraints, setCostFunctionForClosestCP,
setCostFunction, setExplicitAEs, setGeneralConstraints, setLinearConstraints,
setObjfunTemplatePath, setParametersForAugSys, setParameters, setVariables, 
setVarsForAugSysWithoutRange, setVarsForClosestPontWithoutRange, 
setVarsWithoutRange; addConstraint := proc (NewConstraint::{0 = term, 0 < term
}) local OldConstraints; OldConstraints := NLPinModule["Constraints"]; if not
type(OldConstraints,{list({0 = term, 0 < term}), []}) then error "\"Constraint\
s\" does not exist (before adding any constaint, \"Constraints\" must exist ev\
en if it is empty) or was defining wrong" end if; NLPinModule["Constraints"] 
:= [op(OldConstraints), NewConstraint]; return end proc; addConstraints := 
proc (listOfConstraints::{list({0 = term, 0 < term}), []}) local 
NumberOfConstraints, i1; NumberOfConstraints := nops(listOfConstraints); for 
i1 to NumberOfConstraints do addConstraint(listOfConstraints[i1]) end do; 
return end proc; addCriticalPointConstraints := proc () local 
CriticalPointConstraints; CriticalPointConstraints := AugSys2:-
ConstraintInitializer:-CriticalPoint:-getConstraints(); addConstraints(
CriticalPointConstraints); return end proc; addDistanceConstraints := proc (
NewNamesOfVars::list(name), GivenDistance::EvalsToFloat, Sign::EvalsToFloat) 
local NewConstraints, RhsNewConstraint, NumOfNVParams, i, NormR, distSign; if
isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
NumOfNVParams := nops(listOfNVParamsInModule); if nops(NewNamesOfVars) <> 
NumOfNVParams then error "Names of Optimized Variables have to be the same ord\
er as list of NV Parameters in augmented system", nops(listOfNVParamsInModule)
end if; NewConstraints := []; if Sign = -1 then distSign := -1 else distSign 
:= 1 end if; NormR := 0; for i to NumOfNVParams do NormR := NormR+r[i]^2 end 
do; NormR := NormR^(1/2); for i to NumOfNVParams do RhsNewConstraint := 
listOfNVParamsInModule[i]-NewNamesOfVars[i]+distSign*dist*r[i]/NormR; 
NewConstraints := [op(NewConstraints), 0 = RhsNewConstraint] end do; 
NewConstraints := [op(NewConstraints), 0 < dist-abs(GivenDistance)]; 
addConstraints(NewConstraints); return end proc; addFixPointConstraints := 
proc () local FixPointConstraints; FixPointConstraints := AugSys2:-
ConstraintInitializer:-FixPoint:-getConstraints(); addConstraints(
FixPointConstraints); return end proc; addGeneralConstraints := proc (
listOfConstraints::{list({0 = term, 0 < term}), []}) local NumberOfConstraints
, Rhs, i1; if not type(NLPinModule["Variables"],list(name = EvalsToFloat .. 
EvalsToFloat)) then error "first define list of \"Variables\" in NLP" end if;
if not type(NLPinModule["ExplicitAEs"],list(name = term)) then error 
"first define list of \"ExplicitAEs\" in NLP" end if; NumberOfConstraints := 
nops(listOfConstraints); for i1 to NumberOfConstraints do Rhs := Aux:-
ListOperations:-subsEqnListIntoEqn(NLPinModule["ExplicitAEs"],
listOfConstraints[i1]); Rhs := rhs(Rhs); if Aux:-LinearEqns:-
checkLinearityOfIn(Rhs,map(lhs,NLPinModule["Variables"])) then 
addLinearConstraint(listOfConstraints[i1]) else addConstraint(
listOfConstraints[i1]) end if end do; return end proc; addLinearConstraint :=
proc (NewConstraint::{0 = term, 0 < term}) local OldConstraints; 
OldConstraints := NLPinModule["LinearConstraints"]; if not type(OldConstraints
,{list({0 = term, 0 < term}), []}) then error "\"LinearConstraints\" does not \
exist (before adding any constaint, \"LinearConstraints\" must exist even if i\
t is empty) or was defining wrong" end if; NLPinModule["LinearConstraints"] :=
[op(OldConstraints), NewConstraint]; return end proc; 
addNVConstraintsForTighterApprBox := proc (namesForNewRobustnessBallVars::list
(name), namesForNVToNewRobustnessBallVars::list(name), 
namesForCriticalBaundVars::list(name), namesForNVVars::list(name), 
namesForOptPointVars::list(name), powerOfNewRobustnessBall::integer, 
minDistance::EvalsToFloat, signOfDirectionNV::integer := 1) local 
numberOfNVVars, i, signNV, newRobBallEqn, squareOfNormNV, 
squareOfNormNVToNewRobBall; if isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
numberOfNVVars := nops(listOfNVParamsInModule); if nops(
namesForCriticalBaundVars) <> numberOfNVVars then error 
"list of names for critical manifold variables has to have dimension", 
numberOfNVVars end if; if nops(namesForNewRobustnessBallVars) <> 
numberOfNVVars then error "list of names for robustness ball variables has not\
 the same dimension as list of variables for critical manifold", 
numberOfNVVars end if; if nops(namesForOptPointVars) <> numberOfNVVars then 
error "list of names for optimal point varibles has not the same dimension as \
list of variables for critical manifold", numberOfNVVars end if; if nops(
namesForNVToNewRobustnessBallVars) <> numberOfNVVars then error "list of names\
 for normal vector varibles to new robustness ball has not the same dimension \
as list of variables for critical manifold", numberOfNVVars end if; if nops(
namesForNVVars) <> numberOfNVVars then error "list of names for normal vector \
varibles has not the same dimension as list of variables for critical manifold\
", numberOfNVVars end if; if convert(namesForCriticalBaundVars,set) <> convert
(listOfNVParamsInModule,set) then error 
"list of names for critical manifold variables has to consists of", 
listOfNVParamsInModule end if; if signOfDirectionNV = -1 then signNV := -1 
else signNV := 1 end if; AugSys2:-ConstraintInitializer:-NLP:-
addDistanceConstraints(namesForNewRobustnessBallVars,minDistance,signNV); for
i to numberOfNVVars do AugSys2:-ConstraintInitializer:-NLP:-addConstraint(0 =
namesForNVToNewRobustnessBallVars[i]-2*powerOfNewRobustnessBall*(
namesForNewRobustnessBallVars[i]-namesForOptPointVars[i])^(2*
powerOfNewRobustnessBall-1)) end do; newRobBallEqn := 2; for i to 
numberOfNVVars do newRobBallEqn := newRobBallEqn-(
namesForNewRobustnessBallVars[i]-namesForOptPointVars[i])^(2*
powerOfNewRobustnessBall) end do; AugSys2:-ConstraintInitializer:-NLP:-
addConstraint(0 = newRobBallEqn); squareOfNormNV := 0; 
squareOfNormNVToNewRobBall := 0; for i to numberOfNVVars do 
squareOfNormNVToNewRobBall := squareOfNormNVToNewRobBall+
namesForNVToNewRobustnessBallVars[i]^2; squareOfNormNV := squareOfNormNV+
namesForNVVars[i]^2 end do; for i to numberOfNVVars do AugSys2:-
ConstraintInitializer:-NLP:-addConstraint(0 = 
namesForNVToNewRobustnessBallVars[i]/squareOfNormNVToNewRobBall^(1/2)+
namesForNVVars[i]/squareOfNormNV^(1/2)) end do; return end proc; 
addNVConstraints := proc () local NVConstraints; NVConstraints := AugSys2:-
ConstraintInitializer:-CriticalPoint:-getNVConstraints(); addConstraints(
NVConstraints); return end proc; addVariable := proc (NewVar::(name = 
EvalsToFloat .. EvalsToFloat)) local OldVars; OldVars := NLPinModule[
"Variables"]; if not type(OldVars,{list(name = EvalsToFloat .. EvalsToFloat),
[]}) then error "\"Variables\" does not exist (before adding any Variable, \"V\
ariables\" must exist even if it is empty) or it was defining wrong" end if; 
NLPinModule["Variables"] := [op(OldVars), NewVar]; return end proc; 
addVariableWithoutRange := proc (NewVar::name) addVariable(NewVar = -infinity
.. infinity); return end proc; copyCriticalPointParamsToNLP := proc () local 
CriticalPointParams; CriticalPointParams := AugSys2:-ConstraintInitializer:-
CriticalPoint:-getParameters(); if not type(CriticalPointParams,list(name = 
EvalsToFloat)) then error "first define Parameters in CriticalPoint" end if; 
NLPinModule["Parameters"] := CriticalPointParams; return end proc; 
createNewNLP := proc () NLPinModule := Aux:-NLP:-newNLP(); return end proc; 
getConstraints := proc () return eval(NLPinModule["Constraints"]) end proc; 
getCostFunction := proc () return eval(NLPinModule["CostFunction"]) end proc;
getExplicitAEs := proc () return eval(NLPinModule["ExplicitAEs"]) end proc; 
getLinearConstraints := proc () return eval(NLPinModule["LinearConstraints"])
end proc; getNLP := proc () return eval(NLPinModule) end proc; 
getObjfunTemplatePath := proc () return NLPinModule["ObjfunTemplatePath"] end
proc; getParameters := proc () return eval(NLPinModule["Parameters"]) end proc
; getVariables := proc () return eval(NLPinModule["Variables"]) end proc; 
setBinaryVariables := proc (newBinVars::list(name)) NLPinModule[
"BinaryVariables"] := newBinVars; return end proc; setConstraints := proc (
listOfConstraints::{list({0 = term, 0 < term}), []}) NLPinModule["Constraints"
] := listOfConstraints; return end proc; setCostFunctionForClosestCP := proc (
pointForFindingDistance::list(name = EvalsToFloat)) local costFunction, i, 
newParsNames, Missing, Obsolete; if isAugSysCreatedInModule <> true then error
"first define augmented system with the procedure createAugSys" end if; 
newParsNames := map(lhs,pointForFindingDistance); Missing, Obsolete := Aux:-
ListOperations:-getMissingAndObsoleteNames(newParsNames,listOfNVParamsInModule
); if not (Missing = {} and Obsolete = {} and nops(pointForFindingDistance) =
nops(listOfNVParamsInModule)) then error "point for which we loking the minima\
l distance to the critical boundary must be given in the form [name1=...,name2\
=..., and so on] where names are equal", listOfNVParamsInModule end if; 
costFunction := 0; for i to nops(pointForFindingDistance) do costFunction := (
lhs(pointForFindingDistance[i])-rhs(pointForFindingDistance[i]))^2+
costFunction end do; NLPinModule["CostFunction"] := [eval(costFunction)]; 
return end proc; setCostFunction := proc (costFunction::list(term)) 
NLPinModule["CostFunction"] := costFunction; return end proc; setExplicitAEs 
:= proc (NewExplicitAEs::list(name = term)) NLPinModule["ExplicitAEs"] := 
NewExplicitAEs; return end proc; setGeneralConstraints := proc (
listOfConstraints::{list({0 = term, 0 < term}), []}) NLPinModule[
"LinearConstraints"] := []; NLPinModule["Constraints"] := []; 
addGeneralConstraints(listOfConstraints); return end proc; 
setLinearConstraints := proc (listOfConstraints::{list({0 = term, 0 < term}),
[]}) NLPinModule["LinearConstraints"] := listOfConstraints; return end proc; 
setObjfunTemplatePath := proc (path::string) NLPinModule["ObjfunTemplatePath"]
:= path; return end proc; setParametersForAugSys := proc (NewPars::list(name =
EvalsToFloat)) local NewParsNames, ParNames, Pars, Missing, Obsolete; if 
isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
NewParsNames := map(lhs,NewPars); NewParsNames := [op(NewParsNames), op(
listOfNVParamsInModule)]; Pars := augSysInModule["Parameters"]; ParNames := 
map(lhs,Pars); Missing, Obsolete := Aux:-ListOperations:-
getMissingAndObsoleteNames(NewParsNames,ParNames); if not Missing = {} then 
error "assignments are missing for %1", Missing end if; NLPinModule[
"Parameters"] := NewPars; return end proc; setParameters := proc (NewPars::
list(name = EvalsToFloat)) NLPinModule["Parameters"] := NewPars; return end 
proc; setVariables := proc (NewVars::{list(name = EvalsToFloat .. EvalsToFloat
)}) NLPinModule["Variables"] := NewVars; return end proc; 
setVarsForAugSysWithoutRange := proc () local listOfVars; if 
isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
listOfVars := [op(listOfNVParamsInModule), op(augSysInModule["Variables"])]; 
setVarsWithoutRange(listOfVars); return end proc; 
setVarsForClosestPontWithoutRange := proc () local listOfVars, 
listOfVarsFromAugSys, augSysVars, i1, numOfFixSystemEqns; if 
isAugSysCreatedInModule <> true then error 
"first define augmented system with the procedure createAugSys" end if; 
numOfFixSystemEqns := nops(listOfDynVarsInModule); augSysVars := 
augSysInModule["Variables"]; if isNSType = true then listOfVarsFromAugSys := [
seq(augSysVars[i1],i1 = 1 .. 3*numOfFixSystemEqns+1)] end if; if 
isFlipOrFoldType = true then listOfVarsFromAugSys := [seq(augSysVars[i1],i1 =
1 .. 2*numOfFixSystemEqns)] end if; listOfVars := [op(listOfNVParamsInModule),
op(listOfVarsFromAugSys)]; setVarsWithoutRange(listOfVars); return end proc; 
setVarsWithoutRange := proc (VarsNames::{list(name)}) local NewVars; NewVars 
:= [seq(VarsNames[i1] = -infinity .. infinity,i1 = 1 .. nops(VarsNames))]; 
NLPinModule["Variables"] := NewVars; return end proc; end module; 
runOptForListStartPoints := proc (NPSOLproc, ListOfStartingPoints, 
NumberParamToPrint::integer) local StrartingCriticalPoint, i, imform1, obj1, 
optVec1, numberStartPoints, bestObj, bestImform, bestNumber, bestOptVect, 
bestList, i2; numberStartPoints := nops(ListOfStartingPoints); bestObj := 
infinity; for i to numberStartPoints do AugSys2:-ConstraintInitializer:-
CriticalPoint:-setStartingCriticalPoint(ListOfStartingPoints[i]); 
StrartingCriticalPoint := AugSys2:-ConstraintInitializer:-CriticalPoint:-
getStartingCriticalPoint(); NPSOLproc:-setXVEC(StrartingCriticalPoint); try 
NPSOLproc:-runOpt(); imform1 := NPSOLproc:-getInform(); obj1 := NPSOLproc:-
getObjf(); optVec1 := convert(NPSOLproc:-getXVEC(),list); printf(
"i=%d, inform=%d, obj=%f, opt=[",i,imform1,obj1); for i2 to NumberParamToPrint
-1 do printf("%f, ",optVec1[i2]) end do; if obj1 < bestObj then bestObj := 
obj1; bestImform := imform1; bestNumber := i; bestOptVect := optVec1 end if; 
printf("%f",optVec1[NumberParamToPrint]); printf("], bestNo=%d\n",bestNumber)
catch: printf("i=%d, time limit error\n",i) end try end do; printf("\n"); 
bestList := [bestNumber, bestImform, bestObj, bestOptVect]; return bestList 
end proc; end module; Discretization := module () export getGaussSaidelDTASys;
getGaussSaidelDTASys := proc (aDAESys::DAESys, listOfFuncParts::listlist) 
local aDAESysEqns, newDTASys, numDAESysEqns, i, gPart, hPart, i1, 
aDESysWithSubs, newDynEqns, newDTASysDynVars, newDTASysAlgVars, 
newDTASysAlgEqns, newDTASysExplicitAlgEqns, newDTASysParameters; 
aDESysWithSubs := Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(aDAESys); 
aDAESysEqns := aDESysWithSubs["ODEs"]; numDAESysEqns := nops(aDAESysEqns); if
nops(listOfFuncParts) <> numDAESysEqns then error 
"number of parts of equations in input has to be the same as ODEs" end if; for
i to numDAESysEqns do if nops(listOfFuncParts[i]) <> 2 then error 
"number of parts of equations in every input has to be equal two" end if; 
gPart[i] := listOfFuncParts[i][1]; hPart[i] := listOfFuncParts[i][2] end do; 
for i1 to numDAESysEqns do if 0 <> simplify(rhs(aDAESysEqns[i1])-gPart[i1]+
aDESysWithSubs["DynVars"][i1]*hPart[i1]) then error "parts of equation", i1, 
"defined wrong", rhs(aDAESysEqns[i1]), not equal, gPart[i1]-aDESysWithSubs[
"DynVars"][i1]*hPart[i1] end if end do; newDTASys := table([]); 
newDTASysDynVars := aDESysWithSubs["DynVars"]; newDTASysAlgVars := 
aDESysWithSubs["AlgVars"]; newDTASysAlgEqns := aDESysWithSubs["AEs"]; 
newDTASysExplicitAlgEqns := aDESysWithSubs["ExplicitAEs"]; newDTASysParameters
:= [op(aDESysWithSubs["Parameters"]), TimeStep = 1]; newDynEqns := []; for i 
to nops(newDTASysDynVars) do newDynEqns := [op(newDynEqns), newDTASysDynVars[i
] = subs(newDynEqns,(newDTASysDynVars[i]+TimeStep*gPart[i])/(1+TimeStep*hPart[
i]))] end do; return table(["DynVars" = newDTASysDynVars, "AlgVars" = 
newDTASysAlgVars, "AlgEqns" = newDTASysAlgEqns, "ExplicitAlgEqns" = 
newDTASysExplicitAlgEqns, "Parameters" = newDTASysParameters, "DynEqns" = 
newDynEqns]) end proc; end module; FlipAndFoldNV := module () export 
calculateNV, CreateFlipNVSys, CreateFlipOrFoldNVSys, 
CreateFlipOrFoldWithAlgEqnsNVSys, CreateFlipWithAlgEqnsNVSys, CreateFoldNVSys,
CreateFoldWithAlgEqnsNVSys, getSysOfEqnsForNV; calculateNV := proc (NVSys::
table, FlipOrFlopEigenVector::Vector, ListOfNVparams::list, lineOfData::list(
name = EvalsToFloat)) local NVSysWithParam, solutionsNV, R, normOfR; 
NVSysWithParam := AugSys2:-FlipAndFoldNV:-getSysOfEqnsForNV(NVSys,
FlipOrFlopEigenVector,lineOfData); solutionsNV := fsolve(NVSysWithParam); R :=
Vector(nops(ListOfNVparams),symbol = r); R := subs(solutionsNV,R); normOfR :=
VectorCalculus[Norm](R); if normOfR <> 0 then R := R/VectorCalculus[Norm](R) 
end if; return R end proc; CreateFlipNVSys := proc (aSys::table, 
ListOfNVparams::list(name)) local flipNVSys, radius; if 2 < nargs then radius
:= args[3]; if not (0 <= radius and radius <= 1) then error "Radius of circle \
where eigenvalues have to lie (third input) have to be between 0 to 1" end if
else radius := 1 end if; flipNVSys := CreateFlipOrFoldNVSys(aSys,
ListOfNVparams,-1,radius); return flipNVSys end proc; CreateFlipOrFoldNVSys :=
proc (aSys::DTASys, ListOfNVparams::list(name), inputEigVal::numeric) local 
ListofNames, item, ComplexRightEigSys, ExtendedSysToBeSubs, 
ExtendedSystemEquations, ExtendedSystemVariables, ExtendedSystemParameters, 
DESys, NumOfDynEqns, f_x, f_xTransp, f_xx, f_xalpha, f_p, f_xpTransp, V, U, R,
NewVariables, NewEquations, f_alpha, f_alphaTransp, W, fxx_w, v_fxx_w, 
fxalpha_w, v_fxalpha_w, NumNVparamsToBeSubs, ParsOfModel, VarsOfModel, result,
NumOfDynEqnsInModule, radius; if 3 < nargs then radius := args[4]; if not (0 
<= radius and radius <= 1) then error "Radius of circle where eigenvalues have\
 to lie (forth input) have to be between 0 to 1" end if else radius := 1 end 
if; ListofNames := map(lhs,aSys["Parameters"]); for item in ListOfNVparams 
while true do if not member(item,ListofNames) then error 
"requested normal vector parameter %1 does not exist in model", item end if 
end do; DESys := Aux:-SystemClasses:-subsExplicitAlgEqnsIntoDTASys(aSys); 
NumOfDynEqns := nops(DESys["DynEqns"]); VarsOfModel := DESys["DynVars"]; 
ParsOfModel := DESys["Parameters"]; f_x := Aux:-Derivs:-f_x(DESys["DynEqns"],
DESys["DynVars"]); f_xTransp := LinearAlgebra[Transpose](f_x); f_alpha := Aux
:-Derivs:-f_p(DESys["DynEqns"],ListOfNVparams); f_alphaTransp := LinearAlgebra
[Transpose](f_alpha); f_xx := Aux:-Derivs:-f_xx(DESys["DynEqns"],DESys[
"DynVars"]); f_xalpha := Aux:-Derivs:-f_xp(DESys["DynEqns"],DESys["DynVars"],
ListOfNVparams); V := [seq(v[i1],i1 = 1 .. NumOfDynEqns)]; W := [seq(w[i1],i1
= 1 .. NumOfDynEqns)]; U := [seq(u[i1],i1 = 1 .. NumOfDynEqns)]; R := [seq(r[
i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := [op(U), op(R)]; 
NewEquations := [seq(0 = rhs(DESys["DynEqns"][i1])-DESys["DynVars"][i1],i1 = 1
.. NumOfDynEqns)]; NewVariables := DESys["DynVars"]; ExtendedSystemEquations 
:= NewEquations; ExtendedSystemVariables := NewVariables; NewEquations := 
LinearAlgebra[Multiply](f_x,convert(W,Vector))-LinearAlgebra[Multiply](convert
(W,Vector),inputEigVal*radius); NewEquations := [seq(0 = NewEquations[i1],i1 =
1 .. NumOfDynEqns)]; NewVariables := [op(W)]; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [op(
ExtendedSystemVariables), op(NewVariables)]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(W,Vector)),convert(W,Vector))-1; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = NewEquations]; 
NewEquations := LinearAlgebra[Multiply](f_xTransp,convert(V,Vector))-
LinearAlgebra[Multiply](convert(V,Vector),inputEigVal*radius)+LinearAlgebra[
Multiply](convert(W,Vector),g1); NewEquations := [seq(0 = NewEquations[i1],i1
= 1 .. NumOfDynEqns)]; NewVariables := [op(V), g1]; ExtendedSystemEquations :=
[op(ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [
op(ExtendedSystemVariables), op(NewVariables)]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(V,Vector)),convert(W,Vector))-1; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = NewEquations]; 
fxx_w := Aux:-TensProd:-Tijk_xj(f_xx,W); v_fxx_w := Aux:-TensProd:-xi_Aij(V,
fxx_w); NewEquations := LinearAlgebra[Multiply](f_xTransp,convert(U,Vector))-
convert(U,Vector)+convert(v_fxx_w,Vector); NewEquations := [seq(0 = 
NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; NewVariables := U; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
fxalpha_w := Aux:-TensProd:-Tijk_xj(f_xalpha,W); v_fxalpha_w := Aux:-TensProd
:-xi_Tij(V,fxalpha_w); NewEquations := LinearAlgebra[Multiply](f_alphaTransp,
convert(U,Vector))+convert(v_fxalpha_w,Vector)-convert(R,Vector); NewEquations
:= [seq(0 = NewEquations[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables :=
R; ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)];
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
ExtendedSysToBeSubs := table(); ExtendedSysToBeSubs["Equations"] := 
ExtendedSystemEquations; ExtendedSysToBeSubs["Variables"] := 
ExtendedSystemVariables; ExtendedSysToBeSubs["Parameters"] := DESys[
"Parameters"]; NumNVparamsToBeSubs := nops(ListOfNVparams); result := module (
) local ExtendedSysInModule, NumNVparamsInModule, JacInModule, 
VarsOfModelInModule, ParsOfModelInModule; export getAugSys, getSys, getNVsys,
getJac, getFlipOrFoldEigenVector, getFoldEigenVector, getFlipEigenVector; 
ExtendedSysInModule := copy(ExtendedSysToBeSubs); NumOfDynEqnsInModule := 
NumOfDynEqns; NumNVparamsInModule := NumNVparamsToBeSubs; VarsOfModelInModule
:= VarsOfModel; ParsOfModelInModule := ParsOfModel; JacInModule := copy(f_x);
getSys := proc () return eval(ExtendedSysInModule) end proc; getJac := proc ()
return eval(JacInModule) end proc; getNVsys := proc () local NVsys, NVvars; 
NVsys := ExtendedSysInModule["Equations"][2*NumOfDynEqnsInModule+2 .. 4*
NumOfDynEqnsInModule+2+NumNVparamsInModule]; NVvars := ExtendedSysInModule[
"Variables"][2*NumOfDynEqnsInModule+1 .. 4*NumOfDynEqnsInModule+1+
NumNVparamsInModule]; return table(["Equations" = NVsys, "Variables" = NVvars]
) end proc; getFlipOrFoldEigenVector := proc (lineOfDataInit::list(name = 
EvalsToFloat), inputEigVal) local ListOfNamesInlineOfDataInit, 
ListOfVarsAndPars, ListOfUnassigneds, fixPointJac0, eigVectors, i1, i2, 
flipEigenVal, flipEigenVectors, flipEigenVector, paramList; 
ListOfNamesInlineOfDataInit := map(lhs,lineOfDataInit); paramList := map(lhs,
ParsOfModelInModule); ListOfVarsAndPars := [op(VarsOfModelInModule), op(
paramList)]; ListOfUnassigneds := []; for item in ListOfVarsAndPars while true
do if not member(item,ListOfNamesInlineOfDataInit) then ListOfUnassigneds := [
op(ListOfUnassigneds), item] end if end do; if not ListOfUnassigneds = [] then
error "variables and parameters %1 must be assigned by 1st argument", 
ListOfUnassigneds end if; fixPointJac0 := evalf(subs(lineOfDataInit,getJac()))
; eigVectors := LinearAlgebra[Eigenvectors](fixPointJac0,output = ('list')); 
flipEigenVal := 0; for i1 while i1 <= nops(eigVectors) and flipEigenVal = 0 do
if abs(Re(eigVectors[i1,1])-inputEigVal) < 1/10000 and abs(Im(eigVectors[i1,1]
)) < 1/10000 then flipEigenVectors := eigVectors[i1,3]; flipEigenVal := 
inputEigVal end if end do; if flipEigenVal = 0 then error 
"There is no point with eigenvalue=", inputEigVal end if; flipEigenVector := 
flipEigenVectors[1]; for i2 to LinearAlgebra[Dimension](flipEigenVector) do if
1/100000 < abs(Im(flipEigenVector[i2])) then error 
"Eigenvector corresponding to eigenvalue=", inputEigVal, 
"is not a real vector" end if end do; flipEigenVector := map(Re,
flipEigenVector); if LinearAlgebra[VectorNorm](flipEigenVector,2) <> 1 then 
flipEigenVector := flipEigenVector/LinearAlgebra[VectorNorm](flipEigenVector,2
) end if; return flipEigenVector end proc; getFlipEigenVector := proc (
lineOfDataInit::list(name = EvalsToFloat)) local flipEigenVector, radius; if 1
< nargs then radius := args[2]; if not (0 <= radius and radius <= 1) then 
error "Radius of circle where eigenvalues have to lie (second input) have to b\
e between 0 to 1" end if else radius := 1 end if; flipEigenVector := 
getFlipOrFoldEigenVector(lineOfDataInit,-radius); return flipEigenVector end 
proc; getFoldEigenVector := proc (lineOfDataInit::list(name = EvalsToFloat)) 
local foldEigenVector, radius; if 1 < nargs then radius := args[2]; if not (0
<= radius and radius <= 1) then error "Radius of circle where eigenvalues have\
 to lie (second input) have to be between 0 to 1" end if else radius := 1 end
if; foldEigenVector := getFlipOrFoldEigenVector(lineOfDataInit,radius); return
foldEigenVector end proc; getAugSys := proc () local AugSysEqns, AugSysVars; 
AugSysEqns := ExtendedSysInModule["Equations"][1 .. 2*NumOfDynEqnsInModule+1];
AugSysVars := ExtendedSysInModule["Variables"][1 .. 2*NumOfDynEqnsInModule]; 
return table(["Equations" = AugSysEqns, "Variables" = AugSysVars, "Parameters"
= ParsOfModel]) end proc end module; return result end proc; 
CreateFlipOrFoldWithAlgEqnsNVSys := proc (aSys::DTASys, ListOfNVparams::list(
name), inputEigVal::numeric) local ListofNames, item, ComplexRightEigSys, 
ExtendedSysToBeSubs, ExtendedSystemEquations, ExtendedSystemVariables, 
ExtendedSystemParameters, DESys, NumOfEqnsOfModel, f_x, f_xTransp, f_xx, 
f_xalpha, f_p, f_xpTransp, V, U, R, NewVariables, NewEquations, f_alpha, 
f_alphaTransp, W, fxx_w, v_fxx_w, fxalpha_w, v_fxalpha_w, NumNVparamsToBeSubs,
ParsOfModel, VarsOfModel, result, NumOfDynEqnsInModule, radius, AlgEqnsOfModel
, AlgVarsOfModel, DynEqnsOfModel, CMatr, i1; if 3 < nargs then radius := args[
4]; if not (0 <= radius and radius <= 1) then error "Radius of circle where ei\
genvalues have to lie (forth input) have to be between 0 to 1" end if else 
radius := 1 end if; ListofNames := map(lhs,aSys["Parameters"]); for item in 
ListOfNVparams while true do if not member(item,ListofNames) then error 
"requested normal vector parameter %1 does not exist in model", item end if 
end do; DESys := Aux:-SystemClasses:-subsExplicitAlgEqnsIntoDTASys(aSys); 
DynEqnsOfModel := DESys["DynEqns"]; AlgEqnsOfModel := DESys["AlgEqns"]; 
NumOfEqnsOfModel := nops(DynEqnsOfModel)+nops(AlgEqnsOfModel); VarsOfModel :=
DESys["DynVars"]; ParsOfModel := DESys["Parameters"]; AlgVarsOfModel := DESys[
"AlgVars"]; if nops(AlgEqnsOfModel) = 0 then return CreateFlipOrFoldNVSys(aSys
,ListOfNVparams,inputEigVal,radius) end if; f_x := Aux:-Derivs:-f_x([op(
DynEqnsOfModel), op(AlgEqnsOfModel)],[op(VarsOfModel), op(AlgVarsOfModel)]); 
f_xTransp := LinearAlgebra[Transpose](f_x); f_alpha := Aux:-Derivs:-f_p([op(
DynEqnsOfModel), op(AlgEqnsOfModel)],ListOfNVparams); f_alphaTransp := 
LinearAlgebra[Transpose](f_alpha); f_xx := Aux:-Derivs:-f_xx([op(
DynEqnsOfModel), op(AlgEqnsOfModel)],[op(VarsOfModel), op(AlgVarsOfModel)]); 
f_xalpha := Aux:-Derivs:-f_xp([op(DynEqnsOfModel), op(AlgEqnsOfModel)],[op(
VarsOfModel), op(AlgVarsOfModel)],ListOfNVparams); V := [seq(v[i1],i1 = 1 .. 
NumOfEqnsOfModel)]; W := [seq(w[i1],i1 = 1 .. NumOfEqnsOfModel)]; U := [seq(u[
i1],i1 = 1 .. NumOfEqnsOfModel)]; R := [seq(r[i1],i1 = 1 .. nops(
ListOfNVparams))]; NewVariables := [op(U), op(R)]; NewEquations := [seq(0 = 
rhs(DESys["DynEqns"][i1])-DESys["DynVars"][i1],i1 = 1 .. nops(DESys["DynEqns"]
))]; NewVariables := DESys["DynVars"]; ExtendedSystemEquations := NewEquations
; ExtendedSystemVariables := NewVariables; NewEquations := [seq(0 = rhs(DESys[
"AlgEqns"][i1]),i1 = 1 .. nops(DESys["AlgEqns"]))]; NewVariables := DESys[
"AlgVars"]; ExtendedSystemEquations := [op(ExtendedSystemEquations), op(
NewEquations)]; ExtendedSystemVariables := [op(ExtendedSystemVariables), op(
NewVariables)]; CMatr := Matrix(NumOfEqnsOfModel); for i1 to nops(DESys[
"DynEqns"]) do CMatr[i1,i1] := 1 end do; NewEquations := LinearAlgebra[
Multiply](f_x,convert(W,Vector))-LinearAlgebra[Multiply](CMatr,LinearAlgebra[
Multiply](convert(W,Vector),inputEigVal*radius)); NewEquations := [seq(0 = 
NewEquations[i1],i1 = 1 .. NumOfEqnsOfModel)]; NewVariables := [op(W)]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
NewEquations := LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W,
Vector)),convert(W,Vector))-1; ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(V,Vector))-LinearAlgebra[Multiply](CMatr,
LinearAlgebra[Multiply](convert(V,Vector),inputEigVal*radius))+LinearAlgebra[
Multiply](convert(W,Vector),g1); NewEquations := [seq(0 = NewEquations[i1],i1
= 1 .. NumOfEqnsOfModel)]; NewVariables := [op(V), g1]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
NewEquations := LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(V,
Vector)),convert(W,Vector))-1; ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; fxx_w := Aux:-TensProd:-Tijk_xj(
f_xx,W); v_fxx_w := Aux:-TensProd:-xi_Aij(V,fxx_w); NewEquations := 
LinearAlgebra[Multiply](f_xTransp,convert(U,Vector))-LinearAlgebra[Multiply](
CMatr,convert(U,Vector))+convert(v_fxx_w,Vector); NewEquations := [seq(0 = 
NewEquations[i1],i1 = 1 .. NumOfEqnsOfModel)]; NewVariables := U; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
fxalpha_w := Aux:-TensProd:-Tijk_xj(f_xalpha,W); v_fxalpha_w := Aux:-TensProd
:-xi_Tij(V,fxalpha_w); NewEquations := LinearAlgebra[Multiply](f_alphaTransp,
convert(U,Vector))+convert(v_fxalpha_w,Vector)-convert(R,Vector); NewEquations
:= [seq(0 = NewEquations[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables :=
R; ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)];
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
ExtendedSysToBeSubs := table(); ExtendedSysToBeSubs["Equations"] := 
ExtendedSystemEquations; ExtendedSysToBeSubs["Variables"] := 
ExtendedSystemVariables; ExtendedSysToBeSubs["Parameters"] := DESys[
"Parameters"]; NumNVparamsToBeSubs := nops(ListOfNVparams); result := module (
) local ExtendedSysInModule, NumNVparamsInModule, JacInModule, 
VarsOfModelInModule, ParsOfModelInModule, NumOfEqnsInModule, 
AlgVarsOfModelInModule; export getAugSys, getSys, getNVsys, getJac, 
getFlipOrFoldEigenVector, getFoldEigenVector, getFlipEigenVector; 
ExtendedSysInModule := copy(ExtendedSysToBeSubs); NumOfEqnsInModule := 
NumOfEqnsOfModel; NumNVparamsInModule := NumNVparamsToBeSubs; 
VarsOfModelInModule := VarsOfModel; ParsOfModelInModule := ParsOfModel; 
AlgVarsOfModelInModule := AlgVarsOfModel; JacInModule := copy(f_x); getSys :=
proc () return eval(ExtendedSysInModule) end proc; getJac := proc () return 
eval(JacInModule) end proc; getNVsys := proc () local NVsys, NVvars; NVsys :=
ExtendedSysInModule["Equations"][2*NumOfEqnsInModule+2 .. 4*NumOfEqnsInModule+
2+NumNVparamsInModule]; NVvars := ExtendedSysInModule["Variables"][2*
NumOfEqnsInModule+1 .. 4*NumOfEqnsInModule+1+NumNVparamsInModule]; return 
table(["Equations" = NVsys, "Variables" = NVvars]) end proc; 
getFlipOrFoldEigenVector := proc (lineOfDataInit::list(name = EvalsToFloat), 
inputEigVal) local ListOfNamesInlineOfDataInit, ListOfVarsAndPars, 
ListOfUnassigneds, fixPointJac0, eigVectors, i1, i2, flipEigenVal, 
flipEigenVectors, flipEigenVector, paramList, cForGenEigVec; 
ListOfNamesInlineOfDataInit := map(lhs,lineOfDataInit); paramList := map(lhs,
ParsOfModelInModule); ListOfVarsAndPars := [op(VarsOfModelInModule), op(
AlgVarsOfModelInModule), op(paramList)]; ListOfUnassigneds := []; for item in
ListOfVarsAndPars while true do if not member(item,ListOfNamesInlineOfDataInit
) then ListOfUnassigneds := [op(ListOfUnassigneds), item] end if end do; if 
not ListOfUnassigneds = [] then error 
"variables and parameters %1 must be assigned by 1st argument", 
ListOfUnassigneds end if; fixPointJac0 := evalf(subs(lineOfDataInit,getJac()))
; cForGenEigVec := Matrix(NumOfEqnsInModule); for i2 to nops(
VarsOfModelInModule) do cForGenEigVec[i2,i2] := 1 end do; eigVectors := 
LinearAlgebra[Eigenvectors](fixPointJac0,cForGenEigVec,output = ('list')); 
flipEigenVal := 0; for i1 while i1 <= nops(eigVectors) and flipEigenVal = 0 do
if abs(Re(eigVectors[i1,1])-inputEigVal) < 1/10000 and abs(Im(eigVectors[i1,1]
)) < 1/10000 then flipEigenVectors := eigVectors[i1,3]; flipEigenVal := 
inputEigVal end if end do; if flipEigenVal = 0 then error 
"There is no point with eigenvalue=", inputEigVal end if; flipEigenVector := 
flipEigenVectors[1]; for i2 to LinearAlgebra[Dimension](flipEigenVector) do if
1/100000 < abs(Im(flipEigenVector[i2])) then error 
"Eigenvector corresponding to eigenvalue=", inputEigVal, 
"is not a real vector" end if end do; flipEigenVector := map(Re,
flipEigenVector); if LinearAlgebra[VectorNorm](flipEigenVector,2) <> 1 then 
flipEigenVector := flipEigenVector/LinearAlgebra[VectorNorm](flipEigenVector,2
) end if; return flipEigenVector end proc; getFlipEigenVector := proc (
lineOfDataInit::list(name = EvalsToFloat)) local flipEigenVector, radius; if 1
< nargs then radius := args[2]; if not (0 <= radius and radius <= 1) then 
error "Radius of circle where eigenvalues have to lie (second input) have to b\
e between 0 to 1" end if else radius := 1 end if; flipEigenVector := 
getFlipOrFoldEigenVector(lineOfDataInit,-radius); return flipEigenVector end 
proc; getFoldEigenVector := proc (lineOfDataInit::list(name = EvalsToFloat)) 
local foldEigenVector, radius; if 1 < nargs then radius := args[2]; if not (0
<= radius and radius <= 1) then error "Radius of circle where eigenvalues have\
 to lie (second input) have to be between 0 to 1" end if else radius := 1 end
if; foldEigenVector := getFlipOrFoldEigenVector(lineOfDataInit,radius); return
foldEigenVector end proc; getAugSys := proc () local AugSysEqns, AugSysVars; 
AugSysEqns := ExtendedSysInModule["Equations"][1 .. 2*NumOfEqnsInModule+1]; 
AugSysVars := ExtendedSysInModule["Variables"][1 .. 2*NumOfEqnsInModule]; 
return table(["Equations" = AugSysEqns, "Variables" = AugSysVars, "Parameters"
= ParsOfModel]) end proc end module; return result end proc; 
CreateFlipWithAlgEqnsNVSys := proc (aSys::table, ListOfNVparams::list(name)) 
local flipNVSys, radius; if 2 < nargs then radius := args[3]; if not (0 <= 
radius and radius <= 1) then error "Radius of circle where eigenvalues have to\
 lie (third input) have to be between 0 to 1" end if else radius := 1 end if;
flipNVSys := CreateFlipOrFoldWithAlgEqnsNVSys(aSys,ListOfNVparams,-1,radius);
return flipNVSys end proc; CreateFoldNVSys := proc (aSys::table, 
ListOfNVparams::list(name)) local foldNVSys, radius; if 2 < nargs then radius
:= args[3]; if not (0 <= radius and radius <= 1) then error "Radius of circle \
where eigenvalues have to lie (third input) have to be between 0 to 1" end if
else radius := 1 end if; foldNVSys := CreateFlipOrFoldNVSys(aSys,
ListOfNVparams,1,radius); return foldNVSys end proc; 
CreateFoldWithAlgEqnsNVSys := proc (aSys::table, ListOfNVparams::list(name)) 
local foldNVSys, radius; if 2 < nargs then radius := args[3]; if not (0 <= 
radius and radius <= 1) then error "Radius of circle where eigenvalues have to\
 lie (third input) have to be between 0 to 1" end if else radius := 1 end if;
foldNVSys := CreateFlipOrFoldWithAlgEqnsNVSys(aSys,ListOfNVparams,1,radius); 
return foldNVSys end proc; getSysOfEqnsForNV := proc (FlipOrFlopNVSys::table,
FlipOrFlopEigenVector::Vector, lineOfData::list(name = EvalsToFloat)) local 
NVSys2, WParam, i1, listW, NVSys3; NVSys2 := subs(lineOfData,FlipOrFlopNVSys[
"Equations"]); WParam := Vector(LinearAlgebra[Dimension](FlipOrFlopEigenVector
),symbol = w); listW := []; for i1 to LinearAlgebra[Dimension](
FlipOrFlopEigenVector) do listW := [op(listW), WParam[i1] = 
FlipOrFlopEigenVector[i1]] end do; NVSys3 := subs(listW,NVSys2); return eval(
NVSys3) end proc; end module; HopfNV := module () export CreateHopfNVSys; 
CreateHopfNVSys := proc (aSys::DAESys, ListOfNVparams::list(name)) local 
ListofNames, item, ComplexRightEigSys, ExtendedSysToBeSubs, 
ExtendedSystemEquations, ExtendedSystemVariables, ExtendedSystemParameters, 
DESys, NumOfDynEqns, f_x, f_xTransp, f_xx, f_xalpha, f_p, f_xpTransp, V1, V2,
U, R, NewVariables, NewEquations, f_alpha, f_alphaTransp, W1, W2, fxx_w1, 
v1_fxx_w1, fxx_w2, v2_fxx_w2, fxalpha_w1, v1_fxalpha_w1, fxalpha_w2, 
v2_fxalpha_w2, NumNVparamsToBeSubs, ParsOfModel, VarsOfModel, result, 
NumOfDynEqnsInModule, radius; if 2 < nargs then radius := args[3]; if 0 < 
radius then error "The given baund for eigenvalues should be negative" end if
else radius := 0 end if; ListofNames := map(lhs,aSys["Parameters"]); for item
in ListOfNVparams while true do if not member(item,ListofNames) then error 
"requested normal vector parameter %1 does not exist in model", item end if 
end do; DESys := Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(aSys); 
NumOfDynEqns := nops(DESys["ODEs"]); VarsOfModel := DESys["DynVars"]; 
ParsOfModel := DESys["Parameters"]; f_x := Aux:-Derivs:-f_x(DESys["ODEs"],
DESys["DynVars"]); f_xTransp := LinearAlgebra[Transpose](f_x); f_alpha := Aux
:-Derivs:-f_p(DESys["ODEs"],ListOfNVparams); f_alphaTransp := LinearAlgebra[
Transpose](f_alpha); f_xx := Aux:-Derivs:-f_xx(DESys["ODEs"],DESys["DynVars"])
; f_xalpha := Aux:-Derivs:-f_xp(DESys["ODEs"],DESys["DynVars"],ListOfNVparams)
; V1 := [seq(v1[i1],i1 = 1 .. NumOfDynEqns)]; V2 := [seq(v2[i1],i1 = 1 .. 
NumOfDynEqns)]; W1 := [seq(w1[i1],i1 = 1 .. NumOfDynEqns)]; W2 := [seq(w2[i1],
i1 = 1 .. NumOfDynEqns)]; U := [seq(u[i1],i1 = 1 .. NumOfDynEqns)]; R := [seq(
r[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := [op(U), op(R)]; 
NewEquations := [seq(0 = rhs(DESys["ODEs"][i1]),i1 = 1 .. NumOfDynEqns)]; 
NewVariables := DESys["DynVars"]; ExtendedSystemEquations := NewEquations; 
ExtendedSystemVariables := NewVariables; NewEquations := LinearAlgebra[
Multiply](f_x,convert(W1,Vector))-LinearAlgebra[Multiply](radius,convert(W1,
Vector))+LinearAlgebra[Multiply](sigma,convert(W2,Vector)); NewEquations := [
seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; NewVariables := [sigma, op(
W1), op(W2)]; ExtendedSystemEquations := [op(ExtendedSystemEquations), op(
NewEquations)]; ExtendedSystemVariables := [op(ExtendedSystemVariables), op(
NewVariables)]; NewEquations := LinearAlgebra[Multiply](f_x,convert(W2,Vector)
)-LinearAlgebra[Multiply](sigma,convert(W1,Vector))-LinearAlgebra[Multiply](
radius,convert(W2,Vector)); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 
.. NumOfDynEqns)]; ExtendedSystemEquations := [op(ExtendedSystemEquations), op
(NewEquations)]; NewEquations := LinearAlgebra[Multiply](LinearAlgebra[
Transpose](convert(W1,Vector)),convert(W1,Vector))+LinearAlgebra[Multiply](
LinearAlgebra[Transpose](convert(W2,Vector)),convert(W2,Vector))-1; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = NewEquations]; 
NewEquations := LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W1,
Vector)),convert(W2,Vector)); ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(V1,Vector))-LinearAlgebra[Multiply](radius,convert
(V1,Vector))-LinearAlgebra[Multiply](sigma,convert(V2,Vector))+LinearAlgebra[
Multiply](convert(W1,Vector),g1)-LinearAlgebra[Multiply](convert(W2,Vector),g2
); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; 
NewVariables := [op(V1), op(V2), g1, g2]; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [op(
ExtendedSystemVariables), op(NewVariables)]; NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(V2,Vector))+LinearAlgebra[Multiply](sigma,convert(
V1,Vector))-LinearAlgebra[Multiply](radius,convert(V2,Vector))+LinearAlgebra[
Multiply](convert(W2,Vector),g1)+LinearAlgebra[Multiply](convert(W1,Vector),g2
); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
NewEquations := LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(V1,
Vector)),convert(W1,Vector))+LinearAlgebra[Multiply](LinearAlgebra[Transpose](
convert(V2,Vector)),convert(W2,Vector))-1; ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(V1,Vector)),convert(W2,Vector))-
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(V2,Vector)),convert(
W1,Vector)); ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = 
NewEquations]; fxx_w1 := Aux:-TensProd:-Tijk_xj(f_xx,W1); v1_fxx_w1 := Aux:-
TensProd:-xi_Aij(V1,fxx_w1); fxx_w2 := Aux:-TensProd:-Tijk_xj(f_xx,W2); 
v2_fxx_w2 := Aux:-TensProd:-xi_Aij(V2,fxx_w2); NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(U,Vector))+convert(v1_fxx_w1,Vector)+convert(
v2_fxx_w2,Vector); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. 
NumOfDynEqns)]; NewVariables := U; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [op(
ExtendedSystemVariables), op(NewVariables)]; fxalpha_w1 := Aux:-TensProd:-
Tijk_xj(f_xalpha,W1); v1_fxalpha_w1 := Aux:-TensProd:-xi_Tij(V1,fxalpha_w1); 
fxalpha_w2 := Aux:-TensProd:-Tijk_xj(f_xalpha,W2); v2_fxalpha_w2 := Aux:-
TensProd:-xi_Tij(V2,fxalpha_w2); NewEquations := LinearAlgebra[Multiply](
f_alphaTransp,convert(U,Vector))+convert(v1_fxalpha_w1,Vector)+convert(
v2_fxalpha_w2,Vector)-convert(R,Vector); NewEquations := [seq(0 = NewEquations
[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := R; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
ExtendedSysToBeSubs := table(); ExtendedSysToBeSubs["Equations"] := 
ExtendedSystemEquations; ExtendedSysToBeSubs["Variables"] := 
ExtendedSystemVariables; ExtendedSysToBeSubs["Parameters"] := DESys[
"Parameters"]; NumNVparamsToBeSubs := nops(ListOfNVparams); result := module (
) local ExtendedSysInModule, NumNVparamsInModule, JacInModule, 
VarsOfModelInModule, ParsOfModelInModule; export getAugSys, getJac, 
getHopfEigenValueAndEigenVector, getNVsys, getSys; ExtendedSysInModule := copy
(ExtendedSysToBeSubs); NumOfDynEqnsInModule := NumOfDynEqns; 
NumNVparamsInModule := NumNVparamsToBeSubs; VarsOfModelInModule := VarsOfModel
; ParsOfModelInModule := ParsOfModel; JacInModule := copy(f_x); getSys := proc
() return eval(ExtendedSysInModule) end proc; getJac := proc () return eval(
JacInModule) end proc; getNVsys := proc () local NVsys, NVvars; NVsys := 
ExtendedSysInModule["Equations"][3*NumOfDynEqnsInModule+3 .. 6*
NumOfDynEqnsInModule+4+NumNVparamsInModule]; NVvars := ExtendedSysInModule[
"Variables"][3*NumOfDynEqnsInModule+2 .. 6*NumOfDynEqnsInModule+3+
NumNVparamsInModule]; return table(["Equations" = NVsys, "Variables" = NVvars]
) end proc; getHopfEigenValueAndEigenVector := proc (lineOfDataInit::list(name
= EvalsToFloat)) local ListOfNamesInlineOfDataInit, ListOfVarsAndPars, 
ListOfUnassigneds, fixPointJac0, eigVectors, i1, hopfEigenVal, 
hopfEigenVectors, hopfEigenVector, result, reEigenVal, imEigenVal, p, 
paramList, radius; if 1 < nargs then radius := args[2]; if 0 < radius then 
error "The bound on real part of the Hopf eigenvalue should be negative" end 
if else radius := 0 end if; ListOfNamesInlineOfDataInit := map(lhs,
lineOfDataInit); paramList := map(lhs,ParsOfModelInModule); ListOfVarsAndPars
:= [op(VarsOfModelInModule), op(paramList)]; ListOfUnassigneds := []; for item
in ListOfVarsAndPars while true do if not member(item,
ListOfNamesInlineOfDataInit) then ListOfUnassigneds := [op(ListOfUnassigneds),
item] end if end do; if not ListOfUnassigneds = [] then error 
"variables and parameters %1 must be assigned by 1st argument", 
ListOfUnassigneds end if; fixPointJac0 := evalf(subs(lineOfDataInit,getJac()))
; eigVectors := LinearAlgebra[Eigenvectors](fixPointJac0,output = ('list')); 
hopfEigenVal := 0; for i1 while i1 <= nops(eigVectors) and hopfEigenVal = 0 do
if abs(Re(eigVectors[i1,1])-radius) < 1/10000 and 0 < evalf(Im(eigVectors[i1,1
]),6) then hopfEigenVal := eigVectors[i1,1]; hopfEigenVectors := eigVectors[i1
,3] end if end do; if hopfEigenVal = 0 then error 
"There is no point corresponding Hopf bifurcation" end if; hopfEigenVector :=
AugSys2:-NeimarkSackerNV:-normalizeVector(hopfEigenVector); result := [
hopfEigenVal, hopfEigenVector]; return result end proc; getAugSys := proc () 
local AugSysEqns, AugSysVars; AugSysEqns := ExtendedSysInModule["Equations"][1
.. 3*NumOfDynEqnsInModule+2]; AugSysVars := ExtendedSysInModule["Variables"][1
.. 3*NumOfDynEqnsInModule+1]; return table(["Equations" = AugSysEqns, 
"Variables" = AugSysVars, "Parameters" = ParsOfModel]) end proc end module; 
return result end proc; end module; NeimarkSackerNV := module () export 
calculateNV, CreateNeimarkSackerNVSys, CreateNSWithAlgEqnsNVSys, 
getSysOfEqnsForNV, normalizeVector; calculateNV := proc (NVSys::table, NSeig::
list, ListOfNVparams::list, lineOfData::list(name = EvalsToFloat)) local 
NVSysWithParam, solutionsNV, R, normOfR; NVSysWithParam := AugSys2:-
NeimarkSackerNV:-getSysOfEqnsForNV(NVSys,NSeig,lineOfData); solutionsNV := 
fsolve(NVSysWithParam); R := Vector(nops(ListOfNVparams),symbol = r); R := 
subs(solutionsNV,R); normOfR := VectorCalculus[Norm](R); if normOfR <> 0 then
R := R/VectorCalculus[Norm](R) end if; return R end proc; 
CreateNeimarkSackerNVSys := proc (aSys::DTASys, ListOfNVparams::list(name)) 
local ListofNames, item, ComplexRightEigSys, ExtendedSysToBeSubs, 
ExtendedSystemEquations, ExtendedSystemVariables, ExtendedSystemParameters, 
DESys, NumOfDynEqns, f_x, f_xTransp, f_xx, f_xalpha, f_p, f_xpTransp, V1, V2,
U, R, NewVariables, NewEquations, f_alpha, f_alphaTransp, W1, W2, fxx_w1, 
v1_fxx_w1, fxx_w2, v2_fxx_w2, fxalpha_w1, v1_fxalpha_w1, fxalpha_w2, 
v2_fxalpha_w2, NumNVparamsToBeSubs, ParsOfModel, VarsOfModel, result, 
NumOfDynEqnsInModule, radius; if 2 < nargs then radius := args[3]; if not (0 
<= radius and radius <= 1) then error "Radius of circle where eigenvalues have\
 to lie (third input) have to be between 0 to 1" end if else radius := 1 end 
if; ListofNames := map(lhs,aSys["Parameters"]); for item in ListOfNVparams 
while true do if not member(item,ListofNames) then error 
"requested normal vector parameter %1 does not exist in model", item end if 
end do; DESys := Aux:-SystemClasses:-subsExplicitAlgEqnsIntoDTASys(aSys); 
NumOfDynEqns := nops(DESys["DynEqns"]); VarsOfModel := DESys["DynVars"]; 
ParsOfModel := DESys["Parameters"]; f_x := Aux:-Derivs:-f_x(DESys["DynEqns"],
DESys["DynVars"]); f_xTransp := LinearAlgebra[Transpose](f_x); f_alpha := Aux
:-Derivs:-f_p(DESys["DynEqns"],ListOfNVparams); f_alphaTransp := LinearAlgebra
[Transpose](f_alpha); f_xx := Aux:-Derivs:-f_xx(DESys["DynEqns"],DESys[
"DynVars"]); f_xalpha := Aux:-Derivs:-f_xp(DESys["DynEqns"],DESys["DynVars"],
ListOfNVparams); V1 := [seq(v1[i1],i1 = 1 .. NumOfDynEqns)]; V2 := [seq(v2[i1]
,i1 = 1 .. NumOfDynEqns)]; W1 := [seq(w1[i1],i1 = 1 .. NumOfDynEqns)]; W2 := [
seq(w2[i1],i1 = 1 .. NumOfDynEqns)]; U := [seq(u[i1],i1 = 1 .. NumOfDynEqns)];
R := [seq(r[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := [op(U), op(R
)]; NewEquations := [seq(0 = rhs(DESys["DynEqns"][i1])-DESys["DynVars"][i1],i1
= 1 .. NumOfDynEqns)]; NewVariables := DESys["DynVars"]; 
ExtendedSystemEquations := NewEquations; ExtendedSystemVariables := 
NewVariables; NewEquations := LinearAlgebra[Multiply](f_x,convert(W1,Vector))-
LinearAlgebra[Multiply](radius,LinearAlgebra[Multiply](convert(W1,Vector),cos(
p)))+LinearAlgebra[Multiply](radius,LinearAlgebra[Multiply](convert(W2,Vector)
,sin(p))); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)];
NewVariables := [p, op(W1), op(W2)]; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [op(
ExtendedSystemVariables), op(NewVariables)]; NewEquations := LinearAlgebra[
Multiply](f_x,convert(W2,Vector))-LinearAlgebra[Multiply](radius,LinearAlgebra
[Multiply](convert(W1,Vector),sin(p)))-LinearAlgebra[Multiply](radius,
LinearAlgebra[Multiply](convert(W2,Vector),cos(p))); NewEquations := [seq(0 =
NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(W1,Vector)),convert(W1,Vector))+
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W2,Vector)),convert(
W2,Vector))-1; ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = 
NewEquations]; NewEquations := LinearAlgebra[Multiply](LinearAlgebra[Transpose
](convert(W1,Vector)),convert(W2,Vector)); ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(V1,Vector))-LinearAlgebra[Multiply](radius,
LinearAlgebra[Multiply](convert(V1,Vector),cos(p)))-LinearAlgebra[Multiply](
radius,LinearAlgebra[Multiply](convert(V2,Vector),sin(p)))+LinearAlgebra[
Multiply](convert(W1,Vector),g1)-LinearAlgebra[Multiply](convert(W2,Vector),g2
); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; 
NewVariables := [op(V1), op(V2), g1, g2]; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [op(
ExtendedSystemVariables), op(NewVariables)]; NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(V2,Vector))+LinearAlgebra[Multiply](radius,
LinearAlgebra[Multiply](convert(V1,Vector),sin(p)))-LinearAlgebra[Multiply](
radius,LinearAlgebra[Multiply](convert(V2,Vector),cos(p)))+LinearAlgebra[
Multiply](convert(W2,Vector),g1)+LinearAlgebra[Multiply](convert(W1,Vector),g2
); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
NewEquations := sin(p)*LinearAlgebra[Multiply](LinearAlgebra[Transpose](
convert(W1,Vector)),convert(V2,Vector))+cos(p)*LinearAlgebra[Multiply](
LinearAlgebra[Transpose](convert(W2,Vector)),convert(V1,Vector))-cos(p)*
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W1,Vector)),convert(
V2,Vector))+sin(p)*LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W2
,Vector)),convert(V2,Vector)); ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(V1,Vector)),convert(W1,Vector))+
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(V2,Vector)),convert(
W2,Vector))-1; ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = 
NewEquations]; fxx_w1 := Aux:-TensProd:-Tijk_xj(f_xx,W1); v1_fxx_w1 := Aux:-
TensProd:-xi_Aij(V1,fxx_w1); fxx_w2 := Aux:-TensProd:-Tijk_xj(f_xx,W2); 
v2_fxx_w2 := Aux:-TensProd:-xi_Aij(V2,fxx_w2); NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(U,Vector))-convert(U,Vector)+convert(v1_fxx_w1,
Vector)+convert(v2_fxx_w2,Vector); NewEquations := [seq(0 = NewEquations[i1],
i1 = 1 .. NumOfDynEqns)]; NewVariables := U; ExtendedSystemEquations := [op(
ExtendedSystemEquations), op(NewEquations)]; ExtendedSystemVariables := [op(
ExtendedSystemVariables), op(NewVariables)]; fxalpha_w1 := Aux:-TensProd:-
Tijk_xj(f_xalpha,W1); v1_fxalpha_w1 := Aux:-TensProd:-xi_Tij(V1,fxalpha_w1); 
fxalpha_w2 := Aux:-TensProd:-Tijk_xj(f_xalpha,W2); v2_fxalpha_w2 := Aux:-
TensProd:-xi_Tij(V2,fxalpha_w2); NewEquations := LinearAlgebra[Multiply](
f_alphaTransp,convert(U,Vector))+convert(v1_fxalpha_w1,Vector)+convert(
v2_fxalpha_w2,Vector)-convert(R,Vector); NewEquations := [seq(0 = NewEquations
[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := R; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
ExtendedSysToBeSubs := table(); ExtendedSysToBeSubs["Equations"] := 
ExtendedSystemEquations; ExtendedSysToBeSubs["Variables"] := 
ExtendedSystemVariables; ExtendedSysToBeSubs["Parameters"] := DESys[
"Parameters"]; NumNVparamsToBeSubs := nops(ListOfNVparams); result := module (
) local ExtendedSysInModule, NumNVparamsInModule, JacInModule, 
VarsOfModelInModule, ParsOfModelInModule; export getAugSys, getJac, 
getNSEigenValueAndEigenVector, getNVsys, getSys; ExtendedSysInModule := copy(
ExtendedSysToBeSubs); NumOfDynEqnsInModule := NumOfDynEqns; 
NumNVparamsInModule := NumNVparamsToBeSubs; VarsOfModelInModule := VarsOfModel
; ParsOfModelInModule := ParsOfModel; JacInModule := copy(f_x); getSys := proc
() return eval(ExtendedSysInModule) end proc; getJac := proc () return eval(
JacInModule) end proc; getNVsys := proc () local NVsys, NVvars; NVsys := 
ExtendedSysInModule["Equations"][3*NumOfDynEqnsInModule+3 .. 6*
NumOfDynEqnsInModule+4+NumNVparamsInModule]; NVvars := ExtendedSysInModule[
"Variables"][3*NumOfDynEqnsInModule+2 .. 6*NumOfDynEqnsInModule+3+
NumNVparamsInModule]; return table(["Equations" = NVsys, "Variables" = NVvars]
) end proc; getNSEigenValueAndEigenVector := proc (lineOfDataInit::list(name =
EvalsToFloat)) local ListOfNamesInlineOfDataInit, ListOfVarsAndPars, 
ListOfUnassigneds, fixPointJac0, eigVectors, i1, naimarkSackerEigenVal, 
naimarkSackerEigenVectors, naimarkSackerEigenVector, result, reEigenVal, 
imEigenVal, p, paramList, radius; if 1 < nargs then radius := args[2]; if not
(0 <= radius and radius <= 1) then error "Radius of circle where eigenvalues h\
ave to lie (second input) have to be between 0 to 1" end if else radius := 1 
end if; ListOfNamesInlineOfDataInit := map(lhs,lineOfDataInit); paramList := 
map(lhs,ParsOfModelInModule); ListOfVarsAndPars := [op(VarsOfModelInModule), 
op(paramList)]; ListOfUnassigneds := []; for item in ListOfVarsAndPars while 
true do if not member(item,ListOfNamesInlineOfDataInit) then ListOfUnassigneds
:= [op(ListOfUnassigneds), item] end if end do; if not ListOfUnassigneds = []
then error "variables and parameters %1 must be assigned by 1st argument", 
ListOfUnassigneds end if; fixPointJac0 := evalf(subs(lineOfDataInit,getJac()))
; eigVectors := LinearAlgebra[Eigenvectors](fixPointJac0,output = ('list')); 
naimarkSackerEigenVal := 0; for i1 while i1 <= nops(eigVectors) and 
naimarkSackerEigenVal = 0 do if abs(norm(eigVectors[i1,1],2)-radius) < 1/10000
and 0 < evalf(Im(eigVectors[i1,1]),6) then naimarkSackerEigenVal := eigVectors
[i1,1]; naimarkSackerEigenVectors := eigVectors[i1,3] end if end do; if 
naimarkSackerEigenVal = 0 then error 
"There is no point corresponding Neimark Sacker bifurcation" end if; 
naimarkSackerEigenVector := naimarkSackerEigenVectors[1]; reEigenVal := Re(
naimarkSackerEigenVal); imEigenVal := Im(naimarkSackerEigenVal); if reEigenVal
= 0 then p := 1/2*evalf[20](Pi) else if 0 < reEigenVal then p := arctan(
imEigenVal/reEigenVal) else p := arctan(imEigenVal/reEigenVal)+evalf[20](Pi) 
end if end if; naimarkSackerEigenVector := AugSys2:-NeimarkSackerNV:-
normalizeVector(naimarkSackerEigenVector); result := [naimarkSackerEigenVal, p
, naimarkSackerEigenVector]; return result end proc; getAugSys := proc () 
local AugSysEqns, AugSysVars; AugSysEqns := ExtendedSysInModule["Equations"][1
.. 3*NumOfDynEqnsInModule+2]; AugSysVars := ExtendedSysInModule["Variables"][1
.. 3*NumOfDynEqnsInModule+1]; return table(["Equations" = AugSysEqns, 
"Variables" = AugSysVars, "Parameters" = ParsOfModel]) end proc end module; 
return result end proc; CreateNSWithAlgEqnsNVSys := proc (aSys::DTASys, 
ListOfNVparams::list(name)) local ListofNames, item, ComplexRightEigSys, 
ExtendedSysToBeSubs, ExtendedSystemEquations, ExtendedSystemVariables, 
ExtendedSystemParameters, DESys, NumOfEqnsOfModel, f_x, f_xTransp, f_xx, 
f_xalpha, f_p, f_xpTransp, V1, V2, U, R, NewVariables, NewEquations, f_alpha,
f_alphaTransp, W1, W2, fxx_w1, v1_fxx_w1, fxx_w2, v2_fxx_w2, fxalpha_w1, 
v1_fxalpha_w1, fxalpha_w2, v2_fxalpha_w2, NumNVparamsToBeSubs, ParsOfModel, 
VarsOfModel, result, radius, AlgEqnsOfModel, AlgVarsOfModel, DynEqnsOfModel, 
CMatr, i1; if 2 < nargs then radius := args[3]; if not (0 <= radius and radius
<= 1) then error "Radius of circle where eigenvalues have to lie (third input)\
 have to be between 0 to 1" end if else radius := 1 end if; ListofNames := map
(lhs,aSys["Parameters"]); for item in ListOfNVparams while true do if not 
member(item,ListofNames) then error 
"requested normal vector parameter %1 does not exist in model", item end if 
end do; DESys := Aux:-SystemClasses:-subsExplicitAlgEqnsIntoDTASys(aSys); 
DynEqnsOfModel := DESys["DynEqns"]; AlgEqnsOfModel := DESys["AlgEqns"]; 
NumOfEqnsOfModel := nops(DynEqnsOfModel)+nops(AlgEqnsOfModel); VarsOfModel :=
DESys["DynVars"]; ParsOfModel := DESys["Parameters"]; AlgVarsOfModel := DESys[
"AlgVars"]; if nops(AlgEqnsOfModel) = 0 then return CreateNeimarkSackerNVSys(
aSys,ListOfNVparams,radius) end if; f_x := Aux:-Derivs:-f_x([op(DynEqnsOfModel
), op(AlgEqnsOfModel)],[op(VarsOfModel), op(AlgVarsOfModel)]); f_xTransp := 
LinearAlgebra[Transpose](f_x); f_alpha := Aux:-Derivs:-f_p([op(DynEqnsOfModel)
, op(AlgEqnsOfModel)],ListOfNVparams); f_alphaTransp := LinearAlgebra[
Transpose](f_alpha); f_xx := Aux:-Derivs:-f_xx([op(DynEqnsOfModel), op(
AlgEqnsOfModel)],[op(VarsOfModel), op(AlgVarsOfModel)]); f_xalpha := Aux:-
Derivs:-f_xp([op(DynEqnsOfModel), op(AlgEqnsOfModel)],[op(VarsOfModel), op(
AlgVarsOfModel)],ListOfNVparams); V1 := [seq(v1[i1],i1 = 1 .. NumOfEqnsOfModel
)]; V2 := [seq(v2[i1],i1 = 1 .. NumOfEqnsOfModel)]; W1 := [seq(w1[i1],i1 = 1 
.. NumOfEqnsOfModel)]; W2 := [seq(w2[i1],i1 = 1 .. NumOfEqnsOfModel)]; U := [
seq(u[i1],i1 = 1 .. NumOfEqnsOfModel)]; R := [seq(r[i1],i1 = 1 .. nops(
ListOfNVparams))]; NewVariables := [op(U), op(R)]; NewEquations := [seq(0 = 
rhs(DESys["DynEqns"][i1])-DESys["DynVars"][i1],i1 = 1 .. nops(DESys["DynEqns"]
))]; NewVariables := DESys["DynVars"]; ExtendedSystemEquations := NewEquations
; ExtendedSystemVariables := NewVariables; NewEquations := [seq(0 = rhs(DESys[
"AlgEqns"][i1]),i1 = 1 .. nops(DESys["AlgEqns"]))]; NewVariables := DESys[
"AlgVars"]; ExtendedSystemEquations := [op(ExtendedSystemEquations), op(
NewEquations)]; ExtendedSystemVariables := [op(ExtendedSystemVariables), op(
NewVariables)]; CMatr := Matrix(NumOfEqnsOfModel); for i1 to nops(DESys[
"DynEqns"]) do CMatr[i1,i1] := 1 end do; NewEquations := LinearAlgebra[
Multiply](f_x,convert(W1,Vector))-LinearAlgebra[Multiply](CMatr,LinearAlgebra[
Multiply](radius,LinearAlgebra[Multiply](convert(W1,Vector),cos(p))))+
LinearAlgebra[Multiply](CMatr,LinearAlgebra[Multiply](radius,LinearAlgebra[
Multiply](convert(W2,Vector),sin(p)))); NewEquations := [seq(0 = NewEquations[
i1],i1 = 1 .. NumOfEqnsOfModel)]; NewVariables := [p, op(W1), op(W2)]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
NewEquations := LinearAlgebra[Multiply](f_x,convert(W2,Vector))-LinearAlgebra[
Multiply](CMatr,LinearAlgebra[Multiply](radius,LinearAlgebra[Multiply](convert
(W1,Vector),sin(p))))-LinearAlgebra[Multiply](CMatr,LinearAlgebra[Multiply](
radius,LinearAlgebra[Multiply](convert(W2,Vector),cos(p)))); NewEquations := [
seq(0 = NewEquations[i1],i1 = 1 .. NumOfEqnsOfModel)]; ExtendedSystemEquations
:= [op(ExtendedSystemEquations), op(NewEquations)]; NewEquations := 
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W1,Vector)),convert(
W1,Vector))+LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W2,Vector
)),convert(W2,Vector))-1; ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(W1,Vector)),convert(W2,Vector)); 
ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = NewEquations]; 
NewEquations := LinearAlgebra[Multiply](f_xTransp,convert(V1,Vector))-
LinearAlgebra[Multiply](CMatr,LinearAlgebra[Multiply](radius,LinearAlgebra[
Multiply](convert(V1,Vector),cos(p))))-LinearAlgebra[Multiply](CMatr,
LinearAlgebra[Multiply](radius,LinearAlgebra[Multiply](convert(V2,Vector),sin(
p))))+LinearAlgebra[Multiply](convert(W1,Vector),g1)-LinearAlgebra[Multiply](
convert(W2,Vector),g2); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. 
NumOfEqnsOfModel)]; NewVariables := [op(V1), op(V2), g1, g2]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
NewEquations := LinearAlgebra[Multiply](f_xTransp,convert(V2,Vector))+
LinearAlgebra[Multiply](CMatr,LinearAlgebra[Multiply](radius,LinearAlgebra[
Multiply](convert(V1,Vector),sin(p))))-LinearAlgebra[Multiply](CMatr,
LinearAlgebra[Multiply](radius,LinearAlgebra[Multiply](convert(V2,Vector),cos(
p))))+LinearAlgebra[Multiply](convert(W2,Vector),g1)+LinearAlgebra[Multiply](
convert(W1,Vector),g2); NewEquations := [seq(0 = NewEquations[i1],i1 = 1 .. 
NumOfEqnsOfModel)]; ExtendedSystemEquations := [op(ExtendedSystemEquations), 
op(NewEquations)]; NewEquations := sin(p)*LinearAlgebra[Multiply](
LinearAlgebra[Transpose](convert(W1,Vector)),LinearAlgebra[Multiply](CMatr,
convert(V2,Vector)))+cos(p)*LinearAlgebra[Multiply](LinearAlgebra[Transpose](
convert(W2,Vector)),LinearAlgebra[Multiply](CMatr,convert(V1,Vector)))-cos(p)*
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(W1,Vector)),
LinearAlgebra[Multiply](CMatr,convert(V2,Vector)))+sin(p)*LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(W2,Vector)),LinearAlgebra[Multiply]
(CMatr,convert(V2,Vector))); ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](LinearAlgebra[Transpose](convert(V1,Vector)),convert(W1,Vector))+
LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(V2,Vector)),convert(
W2,Vector))-1; ExtendedSystemEquations := [op(ExtendedSystemEquations), 0 = 
NewEquations]; fxx_w1 := Aux:-TensProd:-Tijk_xj(f_xx,W1); v1_fxx_w1 := Aux:-
TensProd:-xi_Aij(V1,fxx_w1); fxx_w2 := Aux:-TensProd:-Tijk_xj(f_xx,W2); 
v2_fxx_w2 := Aux:-TensProd:-xi_Aij(V2,fxx_w2); NewEquations := LinearAlgebra[
Multiply](f_xTransp,convert(U,Vector))-LinearAlgebra[Multiply](CMatr,convert(U
,Vector))+convert(v1_fxx_w1,Vector)+convert(v2_fxx_w2,Vector); NewEquations :=
[seq(0 = NewEquations[i1],i1 = 1 .. NumOfEqnsOfModel)]; NewVariables := U; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
fxalpha_w1 := Aux:-TensProd:-Tijk_xj(f_xalpha,W1); v1_fxalpha_w1 := Aux:-
TensProd:-xi_Tij(V1,fxalpha_w1); fxalpha_w2 := Aux:-TensProd:-Tijk_xj(f_xalpha
,W2); v2_fxalpha_w2 := Aux:-TensProd:-xi_Tij(V2,fxalpha_w2); NewEquations := 
LinearAlgebra[Multiply](f_alphaTransp,convert(U,Vector))+convert(v1_fxalpha_w1
,Vector)+convert(v2_fxalpha_w2,Vector)-convert(R,Vector); NewEquations := [seq
(0 = NewEquations[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := R; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
ExtendedSysToBeSubs := table(); ExtendedSysToBeSubs["Equations"] := 
ExtendedSystemEquations; ExtendedSysToBeSubs["Variables"] := 
ExtendedSystemVariables; ExtendedSysToBeSubs["Parameters"] := DESys[
"Parameters"]; NumNVparamsToBeSubs := nops(ListOfNVparams); result := module (
) local ExtendedSysInModule, NumNVparamsInModule, JacInModule, 
VarsOfModelInModule, ParsOfModelInModule, AlgVarsOfModelInModule, 
NumOfEqnsInModule; export getAugSys, getJac, getNSEigenValueAndEigenVector, 
getNVsys, getSys; ExtendedSysInModule := copy(ExtendedSysToBeSubs); 
NumOfEqnsInModule := NumOfEqnsOfModel; NumNVparamsInModule := 
NumNVparamsToBeSubs; VarsOfModelInModule := VarsOfModel; ParsOfModelInModule 
:= ParsOfModel; AlgVarsOfModelInModule := AlgVarsOfModel; JacInModule := copy(
f_x); getSys := proc () return eval(ExtendedSysInModule) end proc; getJac := 
proc () return eval(JacInModule) end proc; getNVsys := proc () local NVsys, 
NVvars; NVsys := ExtendedSysInModule["Equations"][3*NumOfEqnsInModule+3 .. 6*
NumOfEqnsInModule+4+NumNVparamsInModule]; NVvars := ExtendedSysInModule[
"Variables"][3*NumOfEqnsInModule+2 .. 6*NumOfEqnsInModule+3+
NumNVparamsInModule]; return table(["Equations" = NVsys, "Variables" = NVvars]
) end proc; getNSEigenValueAndEigenVector := proc (lineOfDataInit::list(name =
EvalsToFloat)) local ListOfNamesInlineOfDataInit, ListOfVarsAndPars, 
ListOfUnassigneds, fixPointJac0, eigVectors, i1, naimarkSackerEigenVal, 
naimarkSackerEigenVectors, naimarkSackerEigenVector, result, reEigenVal, 
imEigenVal, p, paramList, radius, cForGenEigVec, i2; if 1 < nargs then radius
:= args[2]; if not (0 <= radius and radius <= 1) then error "Radius of circle \
where eigenvalues have to lie (second input) have to be between 0 to 1" end if
else radius := 1 end if; ListOfNamesInlineOfDataInit := map(lhs,lineOfDataInit
); paramList := map(lhs,ParsOfModelInModule); ListOfVarsAndPars := [op(
VarsOfModelInModule), op(AlgVarsOfModelInModule), op(paramList)]; 
ListOfUnassigneds := []; for item in ListOfVarsAndPars while true do if not 
member(item,ListOfNamesInlineOfDataInit) then ListOfUnassigneds := [op(
ListOfUnassigneds), item] end if end do; if not ListOfUnassigneds = [] then 
error "variables and parameters %1 must be assigned by 1st argument", 
ListOfUnassigneds end if; fixPointJac0 := evalf(subs(lineOfDataInit,getJac()))
; cForGenEigVec := Matrix(NumOfEqnsInModule); for i2 to nops(
VarsOfModelInModule) do cForGenEigVec[i2,i2] := 1 end do; eigVectors := 
LinearAlgebra[Eigenvectors](fixPointJac0,cForGenEigVec,output = ('list')); 
naimarkSackerEigenVal := 0; for i1 while i1 <= nops(eigVectors) and 
naimarkSackerEigenVal = 0 do if abs(norm(eigVectors[i1,1],2)-radius) < 1/10000
and 0 < evalf(Im(eigVectors[i1,1]),6) then naimarkSackerEigenVal := eigVectors
[i1,1]; naimarkSackerEigenVectors := eigVectors[i1,3] end if end do; if 
naimarkSackerEigenVal = 0 then error 
"There is no point corresponding generalized Neimark Sacker bifurcation" end 
if; naimarkSackerEigenVector := naimarkSackerEigenVectors[1]; reEigenVal := Re
(naimarkSackerEigenVal); imEigenVal := Im(naimarkSackerEigenVal); if 
reEigenVal = 0 then p := 1/2*evalf[20](Pi) else if 0 < reEigenVal then p := 
arctan(imEigenVal/reEigenVal) else p := arctan(imEigenVal/reEigenVal)+evalf[20
](Pi) end if end if; naimarkSackerEigenVector := AugSys2:-NeimarkSackerNV:-
normalizeVector(naimarkSackerEigenVector); result := [naimarkSackerEigenVal, p
, naimarkSackerEigenVector]; return result end proc; getAugSys := proc () 
local AugSysEqns, AugSysVars; AugSysEqns := ExtendedSysInModule["Equations"][1
.. 3*NumOfEqnsInModule+2]; AugSysVars := ExtendedSysInModule["Variables"][1 ..
3*NumOfEqnsInModule+1]; return table(["Equations" = AugSysEqns, "Variables" =
AugSysVars, "Parameters" = ParsOfModel]) end proc end module; return result 
end proc; getSysOfEqnsForNV := proc (NVSys::table, NSeig::list, lineOfData::
list(name = EvalsToFloat)) local W1, W2, NVSys2, W1Param, W2Param, i1, i2, 
listW, NVSys3, NVSys4; W1 := map(Re,NSeig[3]); W2 := map(Im,NSeig[3]); NVSys2
:= subs(lineOfData,NVSys["Equations"]); W1Param := Vector(LinearAlgebra[
Dimension](W1),symbol = w1); W2Param := Vector(LinearAlgebra[Dimension](W2),
symbol = w2); listW := []; for i1 to LinearAlgebra[Dimension](W1) do listW :=
[op(listW), W1Param[i1] = W1[i1]] end do; for i2 to LinearAlgebra[Dimension](
W2) do listW := [op(listW), W2Param[i2] = W2[i2]] end do; NVSys3 := subs(listW
,NVSys2); NVSys4 := subs(p = NSeig[2],NVSys3); return eval(NVSys4) end proc; 
normalizeVector := proc (ReqaVec::{Vector(datatype = complex[8]), vector(
complex)}) local NormalizedVec, w1, w2, NewVec, f, x0, aVec; if type(ReqaVec,
vector) then aVec := convert(ReqaVec,Vector) else aVec := copy(ReqaVec) end if
; NormalizedVec := LinearAlgebra[Normalize](aVec,Euclidean); w1 := map(Re,
NormalizedVec); w2 := map(Im,NormalizedVec); f := cos(x)*sin(x)*(LinearAlgebra
[DotProduct](w1,w1)-LinearAlgebra[DotProduct](w2,w2))+(cos(x)^2-sin(x)^2)*
LinearAlgebra[DotProduct](w1,w2); x0 := fsolve(f,x = 0 .. Pi); NewVec := 
LinearAlgebra[ScalarMultiply](NormalizedVec,exp(Complex(1)*x0)); return NewVec
end proc; end module; SaddleNodeNV := module () export CreateSaddleNodeNVSys;
CreateSaddleNodeNVSys := proc (aSys::DAESys, ListOfNVparams::list(name)) local
ListofNames, item, ComplexRightEigSys, ExtendedSysToBeSubs, 
ExtendedSystemEquations, ExtendedSystemVariables, ExtendedSystemParameters, 
DESys, NumOfDynEqns, f_x, f_xTransp, V, R, NewVariables, NewEquations, f_alpha
, f_alphaTransp, radius, NumNVparamsToBeSubs, ParsOfModel, VarsOfModel, result
, NumOfDynEqnsInModule; if 2 < nargs then radius := args[3]; if 0 < radius 
then error "The given baund on eigenvalues should be negative" end if else 
radius := 0 end if; ListofNames := map(lhs,aSys["Parameters"]); for item in 
ListOfNVparams while true do if not member(item,ListofNames) then error 
"requested normal vector parameter %1 does not exist in model", item end if 
end do; DESys := Aux:-SystemClasses:-subsExplicitAEsIntoDAESys(aSys); 
NumOfDynEqns := nops(DESys["ODEs"]); VarsOfModel := DESys["DynVars"]; 
ParsOfModel := DESys["Parameters"]; f_x := Aux:-Derivs:-f_x(DESys["ODEs"],
DESys["DynVars"]); f_xTransp := LinearAlgebra[Transpose](f_x); f_alpha := Aux
:-Derivs:-f_p(DESys["ODEs"],ListOfNVparams); f_alphaTransp := LinearAlgebra[
Transpose](f_alpha); V := [seq(v[i1],i1 = 1 .. NumOfDynEqns)]; R := [seq(r[i1]
,i1 = 1 .. nops(ListOfNVparams))]; NewVariables := [op(R)]; NewEquations := [
seq(0 = rhs(DESys["ODEs"][i1]),i1 = 1 .. NumOfDynEqns)]; NewVariables := DESys
["DynVars"]; ExtendedSystemEquations := NewEquations; ExtendedSystemVariables
:= NewVariables; NewEquations := LinearAlgebra[Multiply](f_xTransp,convert(V,
Vector))-LinearAlgebra[Multiply](radius,convert(V,Vector)); NewEquations := [
seq(0 = NewEquations[i1],i1 = 1 .. NumOfDynEqns)]; NewVariables := [op(V)]; 
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
NewEquations := LinearAlgebra[Multiply](LinearAlgebra[Transpose](convert(V,
Vector)),convert(V,Vector))-1; ExtendedSystemEquations := [op(
ExtendedSystemEquations), 0 = NewEquations]; NewEquations := LinearAlgebra[
Multiply](f_alphaTransp,convert(V,Vector))-convert(R,Vector); NewEquations :=
[seq(0 = NewEquations[i1],i1 = 1 .. nops(ListOfNVparams))]; NewVariables := R;
ExtendedSystemEquations := [op(ExtendedSystemEquations), op(NewEquations)]; 
ExtendedSystemVariables := [op(ExtendedSystemVariables), op(NewVariables)]; 
ExtendedSysToBeSubs := table(); ExtendedSysToBeSubs["Equations"] := 
ExtendedSystemEquations; ExtendedSysToBeSubs["Variables"] := 
ExtendedSystemVariables; ExtendedSysToBeSubs["Parameters"] := DESys[
"Parameters"]; NumNVparamsToBeSubs := nops(ListOfNVparams); result := module (
) local ExtendedSysInModule, NumNVparamsInModule, JacInModule, 
VarsOfModelInModule, ParsOfModelInModule; export getAugSys, getSys, getNVsys,
getJac; ExtendedSysInModule := copy(ExtendedSysToBeSubs); NumOfDynEqnsInModule
:= NumOfDynEqns; NumNVparamsInModule := NumNVparamsToBeSubs; 
VarsOfModelInModule := VarsOfModel; ParsOfModelInModule := ParsOfModel; 
JacInModule := copy(f_x); getSys := proc () return eval(ExtendedSysInModule) 
end proc; getJac := proc () return eval(JacInModule) end proc; getNVsys := 
proc () local NVsys, NVvars; NVsys := ExtendedSysInModule["Equations"][2*
NumOfDynEqnsInModule+2 .. 2*NumOfDynEqnsInModule+1+NumNVparamsInModule]; 
NVvars := ExtendedSysInModule["Variables"][2*NumOfDynEqnsInModule+1 .. 2*
NumOfDynEqnsInModule+NumNVparamsInModule]; return table(["Equations" = NVsys,
"Variables" = NVvars]) end proc; getAugSys := proc () local AugSysEqns, 
AugSysVars; AugSysEqns := ExtendedSysInModule["Equations"][1 .. 2*
NumOfDynEqnsInModule+1]; AugSysVars := ExtendedSysInModule["Variables"][1 .. 2
*NumOfDynEqnsInModule]; return table(["Equations" = AugSysEqns, "Variables" =
AugSysVars, "Parameters" = ParsOfModel]) end proc end module; return result 
end proc; end module; end module;
