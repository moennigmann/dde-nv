%% enter_path_2
% function to enter the path
% you can replace enter_path.m with this function, then you can click
% directly OK in the dialog. You also need to change the call during
% systemInput, this function will give the path as return value, add
% therefor way = enter_path_2.m and remove enter_path.m in systemInput


function y = enter_path_2

f = figure('Visible','off','Position',[360,500,600,120],'Tag','path','Name','Input dialog for maple-path');
% set static head text
headtext = uicontrol('Style','text','Position',[50,40,250,50],'String','Enter the path to maple modules:');
% set input text
editpath = uicontrol('Style','edit','String','/this/is/the/path/to/maple','Position',[50,20,450,25],...
    'Callback',@editMaxReal_Callback);
% set the ok-pushbutton to leave
okpush = uicontrol('Position',[550,20,50,50],'String','OK','Callback',@OKbutton_Callback);
% set visible
f.Visible = 'on';

waitfor(findobj('-regexp','Tag','path'));

function OKbutton_Callback(source,eventdata)
y = editpath.String;     % set input text to outputvariable
close all;                  % close window -> continue
end
end