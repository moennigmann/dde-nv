function [ tau ] = NCSTRdelays( k,~,p )
% gives back delay in 3 CSTR model for use with DDEBIFTool

sumF=p(1);

for i=2:k
    sumF=sumF+p((i-1)*4);
end

tau=5*350/sumF;    
    
end

