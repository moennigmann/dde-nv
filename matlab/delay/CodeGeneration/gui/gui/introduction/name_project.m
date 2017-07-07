%% Introtuction to the gui ask for a name of the project -> to name codegen etc better

function y = intro
f = figure('Visible','off','Position',[360,500,600,170],'Tag','eqn','Name','Input dialog for system dynamics');
% set static head text
uicontrol('Style','text','Position',[50,60,500,60],'String','Welcome to the GUI for your project! Please give your project a name:');
% set input text
editname = uicontrol('Style','edit','String','Project_Name','Position',[50,20,500,25],...
    'Callback',@editName_Callback);
% set the ok-pushbutton to leave
okpush = uicontrol('Position',[550,20,50,50],'String','OK','Callback',@OKbutton_Callback);
f.Visible = 'on';
% y = zeros(1,4);
waitfor(findobj('-regexp','Tag','eqn'));

function editName_Callback(source,eventdata)
    y = editname.String;     % set input text to outputvariable
    close all;                  % close window -> continue
end
end