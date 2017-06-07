% fehler wenn Anzahl verstellt wird wird der wert der bestehenden
% eingebefenster nicht übernommen, müssen also neu eingegeben werden oder
% fenster markieren und enter


function variable_zustand_input
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

%  Create and then hide the UI as it is being constructed.
f = figure('Visible','off','Position',[360,500,450,350],'Tag','Var_z_in');
guidata(f,struct('val',0,'name',[]));
% Construct the components.
hpopup = uicontrol('Style','popupmenu',...
           'String',{1,2,3,4,5,6},...
           'Position',[300,200,100,25],...
           'Callback',@popup_menu_Callback,'Tag','ddmenu','Userdata',struct('val',0));
hx1 = uicontrol('Style','edit','Position',[100,250,100,25],'Callback',@edit1_callback);
hx2 = uicontrol('Style','edit','Position',[100,200,100,25],'Callback',@edit2_callback);
hx3 = uicontrol('Style','edit','Position',[100,150,100,25],'Callback',@edit3_callback);
hx4 = uicontrol('Style','edit','Position',[100,100,100,25],'Callback',@edit4_callback);
hx5 = uicontrol('Style','edit','Position',[100,50,100,25],'Callback',@edit5_callback);
hx6 = uicontrol('Style','edit','Position',[100,0,100,25],'Callback',@edit6_callback);
ht1 = uicontrol('Style','text','Position',[0,250,100,25],'String','x1');
ht2 = uicontrol('Style','text','Position',[0,200,100,25],'String','x2');
ht3 = uicontrol('Style','text','Position',[0,150,100,25],'String','x3');
ht4 = uicontrol('Style','text','Position',[0,100,100,25],'String','x4');
ht5 = uicontrol('Style','text','Position',[0,50,100,25],'String','x5');
ht6 = uicontrol('Style','text','Position',[0,0,100,25],'String','x6');
hstattext = uicontrol('Style','text','Position',[300,250,100,30],'String','Number of states');
hvar = uicontrol('Style','text','Position',[0,300,100,25],'String','Variable');
hname = uicontrol('Style','text','Position',[100,300,100,25],'String','Name');
hpush = uicontrol('Style','pushbutton','String','OK','Position',[300,150,100,25],...
    'Callback',@okbutton_Callback);
% align([hpopup,hpush],'Center','None');

% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hpopup.Units = 'normalized';
hstattext.Units = 'normalized';
hpush.Units = 'normalized';
hx1.Units = 'normalized';
hx2.Units = 'normalized';
hvar.Units = 'normalized';
hname.Units = 'normalized';
% Assign the a name to appear in the window title.
f.Name = 'Input GUI - Give number and names of dynamic variables';

% Move the window to the center of the screen.
movegui(f,'center')

% Make the window visible.
f.Visible = 'on';
% hstattext.Visible = 'on';
hx1.Visible = 'off';
hx2.Visible = 'off';
hx3.Visible = 'off';
hx4.Visible = 'off';
hx5.Visible = 'off';
hx6.Visible = 'off';
ht1.Visible = 'off';
ht2.Visible = 'off';
ht3.Visible = 'off';
ht4.Visible = 'off';
ht5.Visible = 'off';
ht6.Visible = 'off';


function popup_menu_Callback(source,eventdata) 
      % Determine the selected data set.
%       str = get(source, 'String')
      val = get(source,'Value');
%       assignin('base','val',val);
      data.val = val;
      % Set displayed inputs to the selected data set.
      switch val;
      case 1 % User selects 1.
          if strcmp(hx1.Visible,'off')
          hx1.Visible = 'on';
          end
          hx2.Visible = 'off';
          hx3.Visible = 'off';
          hx4.Visible = 'off';
          hx5.Visible = 'off';
          hx6.Visible = 'off';
          if strcmp(hx1.Enable,'off')
          hx1.Enable = 'on';
          end
          hx2.Enable = 'off';
          hx3.Enable = 'off';
          hx4.Enable = 'off';
          hx5.Enable = 'off';
          hx6.Enable = 'off';
          ht1.Visible = 'on';
          ht2.Visible = 'off';
          ht3.Visible = 'off';
          ht4.Visible = 'off';
          ht5.Visible = 'off';
          ht6.Visible = 'off';
      case 2 % User selects 2.
          if strcmp(hx1.Visible,'off')
          hx1.Visible = 'on';
          end
          if strcmp(hx2.Visible,'off')
          hx2.Visible = 'on';
          end
          hx3.Visible = 'off';
          hx4.Visible = 'off';
          hx5.Visible = 'off';
          hx6.Visible = 'off';
          if strcmp(hx1.Enable,'off')
          hx1.Enable = 'on';
          end
          if strcmp(hx2.Enable,'off')
          hx2.Enable = 'on';
          end
          hx3.Enable = 'off';
          hx4.Enable = 'off';
          hx5.Enable = 'off';
          hx6.Enable = 'off';
          ht1.Visible = 'on';
          ht2.Visible = 'on';
          ht3.Visible = 'off';
          ht4.Visible = 'off';
          ht5.Visible = 'off';
          ht6.Visible = 'off';
      case 3 % User selects 3.
          hx1.Visible = 'on';
          hx2.Visible = 'on';
          hx3.Visible = 'on';
          hx4.Visible = 'off';
          hx5.Visible = 'off';
          hx6.Visible = 'off';
          hx1.Enable = 'on';
          hx2.Enable = 'on';
          hx3.Enable = 'on';
          hx4.Enable = 'off';
          hx5.Enable = 'off';
          hx6.Enable = 'off';
          ht1.Visible = 'on';
          ht2.Visible = 'on';
          ht3.Visible = 'on';
          ht4.Visible = 'off';
          ht5.Visible = 'off';
          ht6.Visible = 'off';
      case 4 % User selects 4.
          hx1.Visible = 'on';
          hx2.Visible = 'on';
          hx3.Visible = 'on';
          hx4.Visible = 'on';
          hx5.Visible = 'off';
          hx6.Visible = 'off';
          hx1.Enable = 'on';
          hx2.Enable = 'on';
          hx3.Enable = 'on';
          hx4.Enable = 'on';
          hx5.Enable = 'off';
          hx6.Enable = 'off';
          ht1.Visible = 'on';
          ht2.Visible = 'on';
          ht3.Visible = 'on';
          ht4.Visible = 'on';
          ht5.Visible = 'off';
          ht6.Visible = 'off';
      case 5 % User selects 5.
          hx1.Visible = 'on';
          hx2.Visible = 'on';
          hx3.Visible = 'on';
          hx4.Visible = 'on';
          hx5.Visible = 'on';
          hx6.Visible = 'off';
          hx1.Enable = 'on';
          hx2.Enable = 'on';
          hx3.Enable = 'on';
          hx4.Enable = 'on';
          hx5.Enable = 'on';
          hx6.Enable = 'off';
          ht1.Visible = 'on';
          ht2.Visible = 'on';
          ht3.Visible = 'on';
          ht4.Visible = 'on';
          ht5.Visible = 'on';
          ht6.Visible = 'off';
      case 6 % User selects 6.
          hx1.Visible = 'on';
          hx2.Visible = 'on';
          hx3.Visible = 'on';
          hx4.Visible = 'on';
          hx5.Visible = 'on';
          hx6.Visible = 'on';
          hx1.Enable = 'on';
          hx2.Enable = 'on';
          hx3.Enable = 'on';
          hx4.Enable = 'on';
          hx5.Enable = 'on';
          hx6.Enable = 'on';
          ht1.Visible = 'on';
          ht2.Visible = 'on';
          ht3.Visible = 'on';
          ht4.Visible = 'on';
          ht5.Visible = 'on';
          ht6.Visible = 'on';
      end
      guidata(source,data);
end     
function okbutton_Callback(source,eventdata)
    data = guidata(source);
%     display(data.val);
%     assignin('base','xnum',val);
% data.name = guidata(source);%get(hObject,'String'); % returns contents of edit2 as text
% for i=1:data.val
% display(data.name(i));
% end
assignin('base','names',data.name);
assignin('base','xnum',data.val);
close all;
end

    function edit1_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(1) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edit2_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(2) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edit3_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(3) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edit4_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(4) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edit5_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(5) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edit6_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(6) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
end

