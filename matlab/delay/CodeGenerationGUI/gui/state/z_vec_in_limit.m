%% Script to save names of DynVars to matrix with number and name
% variable xnames containes xi (counted) in first column and name in second
% column

variable_zustand_input;             % opens input dialogue
% uiwait(variable_zustand_input);
waitfor(findobj('-regexp','Tag','Var_z_in'));
xnames{xnum,2}=[];                  % create cell with right size

for i=1:xnum
    xnames(i,:) = [['x' num2str(i)], names(i)]; % write data to cell
end

clear names