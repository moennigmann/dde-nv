%% ask Manifold

function y = askManifold

f = figure('Visible','off','Position',[360,500,350,200],'Tag','Manifold');

% hpushok = uicontrol('Style','pushbutton','String','OK','Position',[350-h,0,150,25],...
%     'Callback',@okbutton_Callback);

headtext = uicontrol('Style','text','Position',[100,160,150,25],'String','Choose the Manifold');
pushFold = uicontrol('Style','checkbox','String','fold','Position',[50,130,80,25],...
    'Callback',@foldButton_Callback);
pushModfold = uicontrol('Style','checkbox','String','modfold','Position',[50,100,80,25],...
    'Callback',@modfoldButton_Callback);
pushHopf = uicontrol('Style','checkbox','String','hopf','Position',[50,70,80,25],...
    'Callback',@hopfButton_Callback);
pushModhopf = uicontrol('Style','checkbox','String','modhopf','Position',[50,40,80,25],...
    'Callback',@modhopfButton_Callback);
pushOK = uicontrol('Style','pushbutton','String','OK','Position',[200,40,80,50],...
    'Callback',@ok_Callback);

f.Visible = 'on';
y = zeros(1,4);
waitfor(findobj('-regexp','Tag','Manifold'));

function foldButton_Callback(source,eventdata)
    if y(1)==0
        y(1) = 1;
    else
        y(1) = 0;
    end
%     close all;
end
function modfoldButton_Callback(source,eventdata)
    if y(2)==0
        y(2) = 1;
    else
        y(2) = 0;
    end
%     close all;
end
function hopfButton_Callback(source,eventdata)
    if y(3)==0
        y(3) = 1;
    else
        y(3) = 0;
    end
%     close all;
end
function modhopfButton_Callback(source,eventdata)
    if y(4)==0
        y(4) = 1;
    else
        y(4) = 0;
    end
%     close all;
end
function ok_Callback(source,eventdata)
    if y==zeros(1,4)
        h = msgbox('Choose at least one','Error','error');
    else
        close all;
    end
end
end