%> @file inputalpha.m
%> @brief opens input window to enter name and value of parameter
%> @return name and value

function varargout = inputalpha(varargin)
% INPUTALPHA MATLAB code for inputalpha.fig
%      INPUTALPHA, by itself, creates a new INPUTALPHA or raises the existing
%      singleton*.
%
%      H = INPUTALPHA returns the handle to a new INPUTALPHA or the handle to
%      the existing singleton*.
%
%      INPUTALPHA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUTALPHA.M with the given input arguments.
%
%      INPUTALPHA('Property','Value',...) creates a new INPUTALPHA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inputalpha_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inputalpha_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inputalpha

% Last Modified by GUIDE v2.5 05-May-2017 11:28:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inputalpha_OpeningFcn, ...
                   'gui_OutputFcn',  @inputalpha_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before inputalpha is made visible.
function inputalpha_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inputalpha (see VARARGIN)
dnum = evalin('base','i');
% Choose default command line output for inputalpha
handles.output = hObject;
handles.data = '_';
handles.val = 0;
% handles.text2 = ['Enter Delay' num2str(2) ':'];
set(handles.text2,'String',['Enter Parameter' num2str(dnum) ' :']);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inputalpha wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inputalpha_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.output = handles.data.del;
% Get default command line output from handles structure
% data = handles.data;
% display('Outputfcn');
% display(handles.data);
% uiresume(handles.figure1);
varargout{1} = handles.data;
varargout{2} = handles.val;
delete(handles.figure1);




function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data = get(hObject,'String'); % returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% guidata(hObject, handles);
handles.output = handles.data;
uiresume(handles.figure1);
% display(handles.data);
% assignin('base','data',handles.data);
% close all;
% delete(handles.figure1);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
  handles.val = str2double(get(hObject,'String')); %returns contents of edit2 as a double
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
