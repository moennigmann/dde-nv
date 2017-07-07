%% enter_path

inputpath;          % opens input dialoge to enter path
waitfor(inputpath);
evalin('base','wayid = fopen(''path.txt'');');  % opens path.txt
evalin('base','way = fscanf(wayid,''%s'');');   % saves path in variable way

clear wayid