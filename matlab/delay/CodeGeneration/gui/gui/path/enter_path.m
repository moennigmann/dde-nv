%> @file enter_path.m
%> @brief code to open input window for entering the path to maple files
%> reads data from path.txt and writes the path directly in the workspace variable way

%% enter_path

inputpath;          % opens input dialoge to enter path
waitfor(inputpath); % wait for inputpath
evalin('base','wayid = fopen(''path.txt'');');  % opens path.txt
evalin('base','way = fscanf(wayid,''%s'');');   % saves path in variable way

clear wayid