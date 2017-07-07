%% Gleichungen einlesen
% to do: 
% anzeige delayedstates funktioniert noch nicht gut für 2 delays oder mehr

% N = xnum;
% N = length(zvec);
% N =2;
% i=1;
fileID = fopen('eqn.txt','w');
% fprintf(fileID,'');
% call for every state the input dialog for the equation
make_samesize;
name_delstates;
eqn = cell(1,xnum);
for i=1:xnum
    xxdot = ['x' num2str(i) 'dot'];
%     inputeqn;           % input dialog - data written to eqn.txt
%     waitfor(inputeqn);
    eqn{i} = enter_eqn;
    waitfor(eqn);
end

% eqnid = fopen('eqn.txt','r');       % open file with data
% 
% for i=1:xnum
% %     eqn(i) = fscanf(eqnid,'%s',[1 i]);
%     eqn{i} = fgetl(eqnid);          % read data from file eqn.txt
% end

for i=1:xnum
    xdot(i,:) = [['x' num2str(i) 'dot'], eqn(i)]; % write data to array
end

clear eqn eqnid fileID xxdot d1