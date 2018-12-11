%> @file eqnabfrage.m
%> @brief calls functions to name delayed states, make inputs the same size
%> for displaying the input in the input dilog, opens input dialog

%% Gleichungen einlesen

% call for every state the input dialog for the equation
name_delstates; % names the delayed states
make_samesize;  % makes entries same size -> needed for enter_eqn
eqn = cell(1,xnum);
for i=1:xnum
    xxdot = ['x' num2str(i) 'dot'];
    eqn{i} = enter_eqn;
    waitfor(eqn);
end

for i=1:xnum
    xdot(i,:) = [['x' num2str(i) 'dot'], eqn(i)]; % write data to array
end

clear eqn eqnid fileID xxdot d1