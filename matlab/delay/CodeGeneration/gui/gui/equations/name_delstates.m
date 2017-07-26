%> @file name_delstates.m
%> @brief creates variable delvars with counted states jointed with delays
%> result looks like : xitauj, i,j: integer

%% Create expressions for delayed states
% delvars = cell(xnum,delnum);

delvars = cell(xnum,delnum);
for j=1:delnum
        for i=1:xnum
            delvars{i,j} = [cell2mat(xnames(i,2)) cell2mat(tauname(j))];
        end
end
