%% enter_equation
% function to ask for the maximal real part for exponential stability at
% modfold and modhopf

function y = enter_eqn

f = figure('Visible','off','Position',[360,500,600,300],'Tag','eqn','Name','Input dialog for system dynamics');
% set static head text
xdotname = evalin('base','xxdot');
xdot = evalin('base','xnames(:,2)');
alvec = evalin('base','alphavec(:,2)');
delvec = evalin('base','del')';
delstates = evalin('base','delvars');
% headtext =
uicontrol('Style','text','Position',[50,200,500,60],'String','Here you are asked to enter the equations for the system. Listed below are the earlier entered states, parameter, delays and delayed states. You can click on them to add them to the equation.');
text = ['Enter the right hand side for ' xdotname ':'];
eqntext = uicontrol('Style','text','Position',[50,30,300,50],'String',text);
% set input text
listheaderstates = uicontrol('Style','text','Position',[50,120,100,60],'String','Possible States');
listheaderparam = uicontrol('Style','text','Position',[200,120,100,60],'String','Possible parameter');
listheaderdelays = uicontrol('Style','text','Position',[350,120,100,60],'String','Possible delays');
listheaderdelstates = uicontrol('Style','text','Position',[500,120,100,60],'String','Possible delayed states');
editeqn = uicontrol('Style','edit','String','equation','Position',[50,20,500,25],...
    'Callback',@editMaxReal_Callback);
% set the ok-pushbutton to leave
okpush = uicontrol('Position',[550,20,50,50],'String','OK','Callback',@OKbutton_Callback);
stateslist = uicontrol('Style','listbox','Position',[50,100,100,50],'String',cell2mat(xdot),'Callback',@stateslist_Callback);
paramlist = uicontrol('Style','listbox','Position',[200,100,100,50],'String',cell2mat(alvec),'Callback',@paramlist_Callback);
delaylist = uicontrol('Style','listbox','Position',[350,100,100,50],'String',cell2mat(delvec),'Callback',@delaylist_Callback);
delstateslist = uicontrol('Style','listbox','Position',[500,100,100,50],'String',cell2mat(delstates),'Callback',@delstateslist_Callback);
% stateslist = uicontrol('Style','listbox','Position',[500,100,100,50],'String',['a'; 'b' ;'c' ;'d'],'Callback',@stateslist_Callback);
f.Visible = 'on';
% y = zeros(1,4);
waitfor(findobj('-regexp','Tag','eqn'));

function stateslist_Callback(source,eventdata)
        editeqn.String = [editeqn.String stateslist.String(stateslist.Value,:)];
end
function paramlist_Callback(source,eventdata)
        editeqn.String = [editeqn.String paramlist.String(paramlist.Value,:)];
end
function delaylist_Callback(source,eventdata)
        editeqn.String = [editeqn.String delaylist.String(delaylist.Value,:)];
end
function delstateslist_Callback(source,eventdata)
        editeqn.String = [editeqn.String delstateslist.String(delstateslist.Value,:)];
end
function editMaxReal_Callback(source,eventdata)
%     y = editMaxReal.String;
end
function OKbutton_Callback(source,eventdata)
%     y = get(handles.editMaxReal,'string');
y = editeqn.String;     % set input text to outputvariable
close all;                  % close window -> continue
end
end