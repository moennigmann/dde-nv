%% Anzahl Totzeiten einlesen

answer = inputdlg('Enter Number of Delaytimes');    % open input dialogue
delnum = str2double(cell2mat(answer));  % get number
tauname(delnum) = {[]};             % create array of right size
for i=1:delnum
    tauname(i) = {['tau' num2str(i)]};  % write data to array
end