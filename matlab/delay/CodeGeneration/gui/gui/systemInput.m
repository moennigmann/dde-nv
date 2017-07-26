%> @file systemInput.m
%> @brief This code opens the input dialogs for defining the system.
%> @param p_name name of the project
%> @param way directory of maple files
%> @param xnames vector with states (names and number)
%> @param xnum number of states
%> @param del contains expressions for delays
%> @param delnum number of delays
%> @param alphavec vector with uncertain parameter (name, velue, number)
%> @param anum number of parameter
%> @param xdot contains equations (rhs) for the states
%> @param choiceMani contains the choice of the manifold (fold,modfold,..)


%% Gui ablauf
% 1.) Initialisierung - Directory eingeben
% 2.) Define System - Variablen und Parameter einlesen
% 3.) Gleichungen einlesen
% Code soll keine beschränkung in anzahl der eingaben haben

%% set matlab
% later anum will contain number of parameters (alpha-vector) and xnum
% number of states

clear
anum = 0;
xnum = 0;

%% Name project
% p_name will contain the name of project -> later use for naming generated
% code

p_name = name_project;

%% Directory einlesen - get directory to maple files
% Pfad in path.txt und WS-variable way speichern
% path saved in variable way

enter_path;

%% Define system
% Dynamische variablen und Namen einlesen
% xnames contains name and number of states
% first column number (xi) second column name
% xnum contains number of states

z_vec_in;

% Anzahl Totzeiten einlesen
% tauname contains counted delays (taui)
% delnum contains number of delays

ask_delaytimes;

% Totzeiten einlesen
% del contains delays in a cell-array

delayabfrage;

% Parameter (alpha-vektor) einlesen
% alphavec contains number, name and value of paramter
% first column number (ai) second name third value
% anum contains number of parameters

a_vec_in;

%% Gleichungen einlesen
% Gleichungen werde zeilenweise in eqn.txt gespeichert und in WS-variable
% xdot
% xdot contains the equations
% first column xidot second column contains rhs of equation

eqnabfrage;

%% Check input
% displays inputs to check and if necessary change them
% pressing in change dialog "ok" will go to case 6 and exit the loop
check = 1;
while check == 1
check_input;
waitfor(findobj('-regexp','Tag','check'));
switch change
    case 1                  % path needs to be changed
        enter_path;
        check = 1;
    case 2                  % change state vector
        z_vec_in_ob;
        check = 1;
    case 3                  % change delaytimes
        ask_delaytimes;
        delayabfrage;
        check = 1;
    case 4                  % change parameter vector
        a_vec_in_ob;
        check = 1;
    case 5                  % change equations
        eqnabfrage;
        check = 1;
    case 6
        break
end
end

clear change check

%% input dialog for additional explicit algebraic equations
% if additional equations are needed (eg:system equations have additional
% parameter)
explicitAEs;


%% ask for kind of mannifold
% choiceMani will be 4-digit array entry 2 and 4 indicate mod

choiceMani = askManifold;

if choiceMani(2) || choiceMani(4)   % if mod fold or modhopf read extra value for max real part
    maxrealpart = ask_maxrealpart;  % open input dialog, save value to maxrealpart
end
