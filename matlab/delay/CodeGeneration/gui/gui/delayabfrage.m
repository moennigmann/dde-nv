%% Input of delays
% array del{} contains delays
% delnum determined by previous function (askdelaytimes)


% del = zeros(delnum,1);
del = cell(1,delnum);    % create cell of right size
for i=1:delnum
    del{i} = inputdelay; % open for every delay input dialogue and save data
%     waitfor(inputdelay);
end


clear names