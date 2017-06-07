% fehler wenn Anzahl verstellt wird wird der wert der bestehenden
% eingebefenster nicht übernommen, müssen also neu eingegeben werden oder
% fenster markieren und enter


function variable_input_2

Nmax = 8;
f = figure('Visible','off','Position',[360,500,450,Nmax*50],'Tag','Var_alpha_in');
guidata(f,struct('val',0,'name',[], 'value',[]));
% Construct the components.
hpopup = uicontrol('Style','popupmenu',...
           'String',{1,2,3,4,5,6,7,8},...
           'Position',[350,250,100,25],...
           'Callback',@popup_menu_Callback,'Tag','ddmenu','Userdata',struct('val',0));   
       
hx(8) = uicontrol('Style','edit','Position',[150,350,100,25],'Callback',@edit8_callback);
hx(7) = uicontrol('Style','edit','Position',[150,300,100,25],'Callback',@edit7_callback);       
hx(6) = uicontrol('Style','edit','Position',[150,250,100,25],'Callback',@edit6_callback);
hx(5) = uicontrol('Style','edit','Position',[150,200,100,25],'Callback',@edit5_callback);
hx(4) = uicontrol('Style','edit','Position',[150,150,100,25],'Callback',@edit4_callback);
hx(3) = uicontrol('Style','edit','Position',[150,100,100,25],'Callback',@edit3_callback);
hx(2) = uicontrol('Style','edit','Position',[150,50,100,25],'Callback',@edit2_callback);
hx(1) = uicontrol('Style','edit','Position',[150,0,100,25],'Callback',@edit1_callback);
hy(8) = uicontrol('Style','edit','Position',[250,350,100,25],'Callback',@edity8_callback);
hy(7) = uicontrol('Style','edit','Position',[250,300,100,25],'Callback',@edity7_callback);       
hy(6) = uicontrol('Style','edit','Position',[250,250,100,25],'Callback',@edity6_callback);
hy(5) = uicontrol('Style','edit','Position',[250,200,100,25],'Callback',@edity5_callback);
hy(4) = uicontrol('Style','edit','Position',[250,150,100,25],'Callback',@edity4_callback);
hy(3) = uicontrol('Style','edit','Position',[250,100,100,25],'Callback',@edity3_callback);
hy(2) = uicontrol('Style','edit','Position',[250,50,100,25],'Callback',@edity2_callback);
hy(1) = uicontrol('Style','edit','Position',[250,0,100,25],'Callback',@edity1_callback);
ht(8) = uicontrol('Style','text','Position',[0,350,150,25],'String','a8 - name - value');
ht(7) = uicontrol('Style','text','Position',[0,300,150,25],'String','a7 - name - value');
ht(6) = uicontrol('Style','text','Position',[0,250,150,25],'String','a6 - name - value');
ht(5) = uicontrol('Style','text','Position',[0,200,150,25],'String','a5 - name - value');
ht(4) = uicontrol('Style','text','Position',[0,150,150,25],'String','a4 - name - value');
ht(3) = uicontrol('Style','text','Position',[0,100,150,25],'String','a3 - name - value');
ht(2) = uicontrol('Style','text','Position',[0,50,150,25],'String','a2 - name - value');
ht(1) = uicontrol('Style','text','Position',[0,0,150,25],'String','a1 - name - value');
hpush = uicontrol('Style','pushbutton','String','OK','Position',[350,150,100,25],...
    'Callback',@okbutton_Callback);
% align([hpopup,hpush],'Center','None');

% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hpopup.Units = 'normalized';
hpush.Units = 'normalized';
% Assign the a name to appear in the window title.
f.Name = 'Input GUI - Give number and names of dynamic variables';

% Move the window to the center of the screen.
movegui(f,'center')

% Make the window visible.
f.Visible = 'on';
for i = 1:Nmax
    hx(i).Visible = 'off';
    ht(i).Visible = 'off';
    hy(i).Visible = 'off';
end


function popup_menu_Callback(source,eventdata) 
      % Determine the selected data set.
%       str = get(source, 'String')
      val = get(source,'Value');
%       assignin('base','val',val);
      data.val = val;
      % Set displayed inputs to the selected data set.
      switch val;
      case 1 % User selects 1.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
      case 2 % User selects 2.
            for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
      case 3 % User selects 3.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
      case 4 % User selects 4.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
      case 5 % User selects 5.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
      case 6 % User selects 6.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
     case 7 % User selects 7.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
     case 8 % User selects 8.
          for i=1:val
              hx(i).Visible = 'on';
              hx(i).Enable = 'on';
              ht(i).Visible = 'on';
              hy(i).Visible = 'on';
              hy(i).Enable = 'on';
          end
          for i=val+1:Nmax
              hx(i).Visible = 'off';
              hx(i).Enable = 'off';
              ht(i).Visible = 'off';
              hy(i).Visible = 'off';
              hy(i).Enable = 'off';
          end
      end
      guidata(source,data);
end     


function okbutton_Callback(source,eventdata)
        data = guidata(source);
%         display(data.val);
%     for i=1:data.val
%         display(data.name(i));
%         display(data.value(i));
%     end
    assignin('base','value',data.value);
    assignin('base','names',data.name);
    assignin('base','anum',data.val);
    close all;
end

    function edit1_callback(source,eventdata)
        data = guidata(source);
        data.name(1) = {source.String};
        guidata(source,data);
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
    function edit7_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(7) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edit8_callback(source,eventdata)
%         data.name1 = get(source,'String');
        data = guidata(source);
        data.name(8) = {source.String};
        guidata(source,data);
%         assignin('base','name',data.name1);
    end
    function edity1_callback(source,eventdata)
        data = guidata(source);
        data.value(1) = str2double(source.String);
        guidata(source,data);
    end
    function edity2_callback(source,eventdata)
        data = guidata(source);
        data.value(2) = str2double(source.String);
        guidata(source,data);
    end
    function edity3_callback(source,eventdata)
        data = guidata(source);
        data.value(3) = str2double(source.String);
        guidata(source,data);
    end
    function edity4_callback(source,eventdata)
        data = guidata(source);
        data.value(4) = str2double(source.String);
        guidata(source,data);
    end
    function edity5_callback(source,eventdata)
        data = guidata(source);
        data.value(5) = str2double(source.String);
        guidata(source,data);
    end
    function edity6_callback(source,eventdata)
        data = guidata(source);
        data.value(6) = str2double(source.String);
        guidata(source,data);
    end
    function edity7_callback(source,eventdata)
        data = guidata(source);
        data.value(7) = str2double(source.String);
        guidata(source,data);
    end
    function edity8_callback(source,eventdata)
        data = guidata(source);
        data.value(8) = str2double(source.String);
        guidata(source,data);
    end
end

