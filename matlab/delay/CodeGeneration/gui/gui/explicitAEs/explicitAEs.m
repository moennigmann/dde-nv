%> @file explicitAEs.m
%> @brief opens input dialog for additionally needed equations
%> if you used in the systems equations parameter wich need a numerical
%> value, you can add theese equations here.
%> @param aes contains the equations

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