%% Create expressions for delayed states
% delvars = cell(xnum,delnum);

delvars = cell(xnum,delnum);
for j=1:delnum
        for i=1:xnum
            delvars{i,j} = [cell2mat(xnames(i,1)) cell2mat(tauname(j))];
        end
end

% if xnum > 10 || delnum >10
% l = 0;
% for i=1:xnum
%     for j=1:delnum
%    d1 = delvars{i,j};
%    if length(d1)>l
%        l = length(d1);
%    end
%     end
% end
% 
% for i=1:xnum
%     for j=1:delnum
%     if length(delvars{i,j})<l
%         for k=1:(l-length(delvars{i,j}))
%         delvars{i,j} = [del{i,j} ' '];
%         end
%     end
%     end
% end 
% end
% dvars = cell(delnum+xnum,1);
% k = 1;
% for i=1:xnum
%     for j=1:delnum
%         dvars{k} = delvars{i,j};
%         k = k+1;
%     end
% end

% clear delvars
% for j=1:delnum
%     for i=1:xnum
%         delvars{i+j-1} = [cell2mat(xnames(i,1)) cell2mat(tauname(j))];
%     end
% end
% delvars{delnum+xnum} = [cell2mat(xnames(xnum,1)) cell2mat(tauname(delnum))];