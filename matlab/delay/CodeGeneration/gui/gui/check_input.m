% display input to check if correct


function check_input
a = evalin('base','xnum');
b = evalin('base','anum');
c = evalin('base','delnum');
% for better eqaution reading
xdote = evalin('base','xdot');
% xdoteqn = cell(a,1);
for i=1:a
xdoteqn{i} = [cell2mat(xdote(i,1)) '=' cell2mat(xdote(i,2))];
end
%
% l = a+b+c;
h = 40;
% h01 = (l+1)*h-35;
% h02 = (l+1)*h-55;
% h01 = (l+1)*h-50;
% h02 = h01-20;
% % h11 = (l+1)*h-85;
% h11 = h02 - 30;
% % h12 = (l+1)*h-(a*20)-85;
% h12 = h11-(a*20);
% % h21 = ((l+1)*h-((a+1)*20))-85;
% h21 = h11-((a+1)*20);
% % h22 = ((l+1)*h-((a+1)*20))-c*20-85;
% h22 = h21 -c*20;
% % h31 = ((l+1)*h-((a+1)*20))-(c+1)*20-85;
% h31 = h21 -(c+1)*20;
% % h32 = ((l+1)*h-((a+1)*20))-(c+1)*20-b*20-85;
% h32 = h31 -b*20;
% % h41 = ((l+1)*h-((a+1)*20))-(c+1)*20-(b+1)*20-85;
% h41 = h31 -(b+1)*20;
% % h42 = ((l+1)*h-((a+1)*20))-(c+1)*20-(b+1)*20-a*20-85;
% h42 = h41 -a*20;
h42 = 20;
h41 = h42 + a*20;
h32 = h41 + 40;
h31 = h32 + b*20;
h22 = h31 + 40;
h21 = h22 + c*20;
h12 = h21 + 40;
h11 = h12 + a*20;
h02 = h11 + 40;
h01 = h02 + 20;
hges = h01 + 30;
%  Create and then hide the UI as it is being constructed.
% f = figure('Visible','off','Position',[360,500,450,350],'Tag','check');
f = figure('Visible','off','Position',[360,500,500,hges],'Tag','check');
% guidata(f,struct('val',0,'name',[]));
% Construct the components.
% hpath1 = uicontrol('Style','text','Position',[0,300,100,25],'String','Path:');
% hpath2 = uicontrol('Style','text','Position',[0,270,100,25],'String',evalin('base','way'));
hpath1 = uicontrol('Style','text','Position',[0,h01,100,25],'String','Path:');
hpath2 = uicontrol('Style','text','Position',[0,h02,300,25],'String',evalin('base','way'));
% hzvec1 = uicontrol('Style','text','Position',[0,240,100,25],'String','States:');
% hzvec2 = uicontrol('Style','text','Position',[30,200,50,40],'String',evalin('base','xnames(:,1)'));
% hzvec3 = uicontrol('Style','text','Position',[70,200,50,40],'String',evalin('base','xnames(:,2)'));
hzvec1 = uicontrol('Style','text','Position',[0,h11,100,25],'String','States:');
hzvec2 = uicontrol('Style','text','Position',[30,h12,50,a*20],'String',evalin('base','xnames(:,1)'));
hzvec3 = uicontrol('Style','text','Position',[80,h12,150,a*20],'String',evalin('base','xnames(:,2)'));
% htauvec1 = uicontrol('Style','text','Position',[0,170,100,25],'String','Delays:');
% htauvec2 = uicontrol('Style','text','Position',[0,140,100,25],'String',evalin('base','del'));
htauvec1 = uicontrol('Style','text','Position',[0,h21,100,25],'String','Delays:');
htauvec2 = uicontrol('Style','text','Position',[0,h22,250,c*20],'String',evalin('base','del'));
% halphavec1 = uicontrol('Style','text','Position',[0,100,100,25],'String','Parameters:');
% halphavec2 = uicontrol('Style','text','Position',[0,70,100,25],'String',evalin('base','alphavec(:,1)'));
% halphavec3 = uicontrol('Style','text','Position',[100,70,50,25],'String',evalin('base','alphavec(:,2)'));
% halphavec4 = uicontrol('Style','text','Position',[150,70,100,25],'String',evalin('base','alphavec(:,3)'));
halphavec1 = uicontrol('Style','text','Position',[0,h31,100,25],'String','Parameters:');
halphavec2 = uicontrol('Style','text','Position',[0,h32,100,b*20],'String',evalin('base','alphavec(:,1)'));
halphavec3 = uicontrol('Style','text','Position',[100,h32,100,b*20],'String',evalin('base','alphavec(:,2)'));
halphavec4 = uicontrol('Style','text','Position',[200,h32,100,b*20],'String',evalin('base','alphavec(:,3)'));
% heqn1 = uicontrol('Style','text','Position',[0,50,100,25],'String','Equations:');
% heqn2 = uicontrol('Style','text','Position',[0,0,50,50],'String',evalin('base','xdot(:,1)'));
% heqn3 = uicontrol('Style','text','Position',[50,0,100,50],'String',evalin('base','xdot(:,2)'));
heqn1 = uicontrol('Style','text','Position',[0,h41,100,25],'String','Equations:');
% heqn2 = uicontrol('Style','text','Position',[0,h42,50,a*20],'String',evalin('base','xdot(:,1)'));
% heqn3 = uicontrol('Style','text','Position',[50,h42,350,a*20],'String',evalin('base','xdot(:,2)'));
heqn2 = uicontrol('Style','listbox','Position',[0,h42,500,a*20],'String',xdoteqn(:));
% heqn3 = uicontrol('Style','text','Position',[50,h42,350,a*20],'String',evalin('base','xdot(:,2)'));
hpushok = uicontrol('Style','pushbutton','String','OK','Position',[350-h,0,150,25],...
    'Callback',@okbutton_Callback);
hpushcp = uicontrol('Style','pushbutton','String','Change Path','Position',[350-h,h01,150,25],...
    'Callback',@cpbutton_Callback);
hpushcs = uicontrol('Style','pushbutton','String','Change States','Position',[350-h,h11,150,25],...
    'Callback',@csbutton_Callback);
hpushcd = uicontrol('Style','pushbutton','String','Change Delays','Position',[350-h,h21,150,25],...
    'Callback',@cdbutton_Callback);
hpushcparam = uicontrol('Style','pushbutton','String','Change Parameters','Position',[350-h,h31,150,25],...
    'Callback',@cparambutton_Callback);
hpushceqn = uicontrol('Style','pushbutton','String','Change Equations','Position',[350-h,h41,150,25],...
    'Callback',@ceqnbutton_Callback);
% align([hpopup,hpush],'Center','None');

% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hstattext.Units = 'normalized';
hpush.Units = 'normalized';
% Assign the a name to appear in the window title.
f.Name = 'Your Inputs to check if correct';

% Move the window to the center of the screen.
movegui(f,'center')

% Make the window visible.
f.Visible = 'on';



function okbutton_Callback(source,eventdata)
    evalin('base','change = 6;');
    close all;
end
function cpbutton_Callback(source,eventdata)
    evalin('base','change = 1;');
    close all;
%     enter_path;
%     check_input;
end
function csbutton_Callback(source,eventdata)
    evalin('base','change = 2;');
    close all;
end
function cdbutton_Callback(source,eventdata)
    evalin('base','change = 3;');
    close all;
end
function cparambutton_Callback(source,eventdata)
    evalin('base','change = 4;');
    close all;
end
function ceqnbutton_Callback(source,eventdata)
    evalin('base','change = 5;');
    close all;
end
end

    

