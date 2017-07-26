%> @file ask_maxrealpart.m
%> @brief opens input dialog for modified manifolds to enter the maximal
%> real part for exponential stability
%> @return maximal real part


%% ask_maxrealpart
% function to ask for the maximal real part for exponential stability at
% modfold and modhopf

function y = ask_maxrealpart

f = figure('Visible','off','Position',[360,500,350,120],'Tag','MaxRealPart');
% set static head text
headtext = uicontrol('Style','text','Position',[50,40,200,50],'String','Enter the maximal real part for exponential stability:');
% set input text
editMaxReal = uicontrol('Style','edit','String','max real part','Position',[50,20,200,25],...
    'Callback',@editMaxReal_Callback);
% set the ok-pushbutton to leave
okpush = uicontrol('Position',[300,20,50,50],'String','OK','Callback',@OKbutton_Callback);


f.Visible = 'on';
% y = zeros(1,4);
waitfor(findobj('-regexp','Tag','MaxRealPart'));

function editMaxReal_Callback(source,eventdata)
%     y = editMaxReal.String;
end
function OKbutton_Callback(source,eventdata)
%     y = get(handles.editMaxReal,'string');
y = editMaxReal.String;     % set input text to outputvariable
close all;                  % close window -> continue
end
end