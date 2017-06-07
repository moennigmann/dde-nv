function varargout = inputx(varargin)
% INPUTX MATLAB code for inputx.fig
%      INPUTX, by itself, creates a new INPUTX or raises the existing
%      singleton*.
%
%      H = INPUTX returns the handle to a new INPUTX or the handle to
%      the existing singleton*.
%
%      INPUTX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUTX.M with the given input arguments.
%
%      INPUTX('Property','Value',...) creates a new INPUTX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inputx_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inputx_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inputx

% Last Modified by GUIDE v2.5 05-May-2017 11:15:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inputx_OpeningFcn, ...
                   'gui_OutputFcn',  @inputx_OutputFcn, ...
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


% --- Executes just before inputx is made visible.
function inputx_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inputx (see VARARGIN)
dnum = evalin('base','i');
% Choose default command line output for inputx
handles.output = hObject;
handles.data = '_';
% handles.text2 = ['Enter Delay' num2str(2) ':'];
set(handles.text2,'String',['Enter State' num2str(dnum) ' :']);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inputx wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inputx_OutputFcn(hObject, eventdata, handles) 
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
% handles.output = handles.data;
uiresume(handles.figure1);
% display(handles.data);
% assignin('base','data',handles.data);
% close all;
% delete(handles.figure1);
