%% Script to save names of Parameters to matrix with number and name

variable_input_2;               % open input dialog
waitfor(findobj('-regexp','Tag','Var_alpha_in'));
alphavec{anum,3}=[];            % create matrix of right size

% write data to martix [ counted parameter, name, value]
for i=1:anum
    alphavec(i,:) = [['alpha' num2str(i)], names(i), value(i)];
end

clear names value answer