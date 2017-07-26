%> @file z_vec_in.m
%> @brief reads in number of states, creats counted states and opens function to read in
%names of states


%% z_vec_in_ob
% einlesen der zustandsvariablen ohne beschränkung der anzahl

answer = inputdlg('Enter Number of States');    % ask for number of states
xnum = str2double(cell2mat(answer));    % save number of states in xnum
xname(xnum) = {[]};             % create array of right size
for i=1:xnum                    % create counted states
    xname(i) = {['x' num2str(i)]};      % xi
end

xnamesp = cell(1,xnum);         % create array of right size
for i=1:xnum                    % array will contain the names of the states
    xnamesp{i} = inputx;        % open input dialog
end

% put the two arrays into one
% column one contains counted states, column two the names
xnames(:,1) = xname(:);
xnames(:,2) = xnamesp(:);

clear xname xnamesp             % delete not needed arrays