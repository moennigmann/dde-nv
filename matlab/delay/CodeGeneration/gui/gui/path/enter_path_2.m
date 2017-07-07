%% enter_path_2
% function to enter the path

function y = enter_path_2

f = figure('Visible','off','Position',[360,500,600,120],'Tag','path','Name','Input dialog for maple-path');
% set static head text
headtext = uicontrol('Style','text','Position',[50,40,250,50],'String','Enter the path to maple modules:');
% set input text
editpath = uicontrol('Style','edit','String','../../../../maple/','Position',[50,20,450,25],...
    'Callback',@editMaxReal_Callback);
% set the ok-pushbutton to leave
okpush = uicontrol('Position',[550,20,50,50],'String','OK','Callback',@OKbutton_Callback);


f.Visible = 'on';
% y = zeros(1,4);
waitfor(findobj('-regexp','Tag','path'));

function editMaxReal_Callback(source,eventdata)
%     y = editMaxReal.String;
end
function OKbutton_Callback(source,eventdata)
%     y = get(handles.editMaxReal,'string');
y = editpath.String;     % set input text to outputvariable
close all;                  % close window -> continue
end
end