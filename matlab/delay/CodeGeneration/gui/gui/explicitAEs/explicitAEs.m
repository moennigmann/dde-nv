%% Input for exlpicit AEs

clear aes
i = 1;
x = 1;
while x
    aes{i} = inputAEs;
    i = i+1;
end
l = length(aes);
aes(l) = [];
l = length(aes);