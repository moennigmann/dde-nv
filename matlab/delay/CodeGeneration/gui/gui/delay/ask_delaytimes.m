%> @file ask_delaytimes.m
%> @brief ask for number of delay times (delnum) and creates counted delays (taui)
%1=1,2,..,delnum

%% Anzahl Totzeiten einlesen

answer = inputdlg('Enter Number of Delaytimes');    % open input dialogue
delnum = str2double(cell2mat(answer));  % get number
tauname(delnum) = {[]};             % create array of right size
for i=1:delnum
    tauname(i) = {['tau' num2str(i)]};  % write data to array
end