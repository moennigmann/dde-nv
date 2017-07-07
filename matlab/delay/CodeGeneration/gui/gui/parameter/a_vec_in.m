%% a_vec_in_ob
% einlesen der zustandsvariablen ohne beschränkung der anzahl
% alphavec wird noch nicht richtig erstellt! probleme bei aval

answer = inputdlg('Enter Number of (uncertain) Parameter');% ask for number of parameter
anum = str2double(cell2mat(answer));        % save number to anum
aname(anum) = {[]};         % create array of right size
for i=1:anum                % create array with counted parameter
    aname(i) = {['alpha' num2str(i)]};      % alphai
end

anamesp = cell(1,anum);         % create arrays of right size
aval = cell(1,anum);
for i=1:anum
    [anamesp{i} aval{i}]= inputalpha;% open input dialog and save to array
end

% put data into one array
alphavec(:,1) = aname(:);
alphavec(:,2) = anamesp(:);
alphavec(:,3) = aval(:);

clear aname anamesp aval            % delete not needed arrays