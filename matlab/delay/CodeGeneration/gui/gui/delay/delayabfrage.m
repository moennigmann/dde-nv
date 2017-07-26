%> @file delayabfrage.m
%> @brief calls function to open inout window and writes data to array del

%% Input of delays
% array del{} contains delays
% delnum determined by previous function (askdelaytimes)

del = cell(1,delnum);    % create cell of right size
for i=1:delnum           % call input dialog delnum times
    del{i} = inputdelay; % open for every delay input dialogue and save data
%     waitfor(inputdelay);
end

clear names