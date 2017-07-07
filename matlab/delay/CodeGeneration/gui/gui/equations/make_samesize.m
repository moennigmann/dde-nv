%% make chars same length

% for xnames 
l = 0;
for i=1:xnum
   d1 = xnames{i,2};
   if length(d1)>l
       l = length(d1);
   end
end

for i=1:xnum
    if length(xnames{i,2})<l
        for j = 1:(l-length(xnames{i,2}))
        xnames{i,2} = [xnames{i,2} ' '];
        end
    end
end

% for delays
l = 0;
for i=1:delnum
   d1 = del{i};
   if length(d1)>l
       l = length(d1);
   end
end

for i=1:delnum
    if length(del{i})<l
        for j = 1:(l-length(del{i}))
        del{i} = [del{i} ' '];
        end
    end
end

% for parameter alpha
l = 0;
for i=1:anum
   d1 = alphavec{i,2};
   if length(d1)>l
       l = length(d1);
   end
end

for i=1:anum
    if length(alphavec{i,2})<l
        for j = 1:(l-length(alphavec{i,2}))
        alphavec{i,2} = [alphavec{i,2} ' '];
        end
    end
end