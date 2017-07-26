%> @file make_samesize.m
%> @brief here the entries will be adjusted by adding spaces to the same
%> size, that is necessary because for displaying them in the input dialog
%> for the equations they need to have the same size (that way cell2mat
%> works, maybe there is a different way that this equalisation is not
%> needed anymore)


%% make chars same length
% needed for enter_eqn otherwise the adding of the entry clicked on eill
% not be added to the eqaution
% for xnames 
l = 0;
% determine the largest length
for i=1:xnum
   d1 = xnames{i,2};
   if length(d1)>l
       l = length(d1);
   end
end
% add spaces to entries shorter than longest
for i=1:xnum
    if length(xnames{i,2})<l
        for j = 1:(l-length(xnames{i,2}))
        xnames{i,2} = [xnames{i,2} ' '];
        end
    end
end

% for delays
% determine longest expression for delays
l = 0;
for i=1:delnum
   d1 = del{i};
   if length(d1)>l
       l = length(d1);
   end
end
% fill shorter expressions with spaces to make them the same length
for i=1:delnum
    if length(del{i})<l
        for j = 1:(l-length(del{i}))
        del{i} = [del{i} ' '];
        end
    end
end

% for parameter alpha
% determine longest entry
l = 0;
for i=1:anum
   d1 = alphavec{i,2};
   if length(d1)>l
       l = length(d1);
   end
end
% make them all the same length
for i=1:anum
    if length(alphavec{i,2})<l
        for j = 1:(l-length(alphavec{i,2}))
        alphavec{i,2} = [alphavec{i,2} ' '];
        end
    end
end


% for delvars
% determine longest entry
l = 0;
for i=1:xnum
    for j=1:delnum
   d1 = delvars{i,j};
   if length(d1)>l
       l = length(d1);
   end
   end
end
% make them all the same length
for i=1:xnum
    for j=1:delnum
    if length(delvars{i,j})<l
        for k = 1:(l-length(delvars{i,j}))
        delvars{i,j} = [delvars{i,j} ' '];
        end
    end
    end
end