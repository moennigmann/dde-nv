%% Gui ablauf
% 1.) Initialisierung - Directory eingeben
% 2.) Define System - Variablen und Parameter einlesen
% 3.) Gleichungen einlesen
% Code ist beschraenkt auf 6 zustaende und 8 parameter

%% set matlab
% later anum will contain number of parameters (alpha-vector) and xnum
% number of states
clear
anum = 0;
xnum = 0;
%% Directory einlesen
% Pfad in path.txt und WS-variable way speichern

enter_path;

%% Define system

% Dynamische variablen und Namen einlesen
% xnames contains name and number of states
% first column number (xi) second column name
% xnum contains number of states

z_vec_in_limit;

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

a_vec_in_limit;

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
l = 0; % used for additional input of explicit AEs
while check == 1
check_input;                            % display inputs for confirmation
waitfor(findobj('-regexp','Tag','check'));
switch change                % switch to input of parameter to change
    case 1          % change path
        enter_path;
        check = 1;
    case 2          % change states
        z_vec_in;
        check = 1;
    case 3          % change delays
        ask_delaytimes;
        delayabfrage;
        check = 1;
    case 4          % change parameter (alpha-vector)
        a_vec_in;
        check = 1;
    case 5          % change equations
        eqnabfrage;
        check = 1;
    case 6          % everything is ok, exit loop
        break
end
end