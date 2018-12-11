%% Function to ask for the input path/directory.
%The path will be saved to the file path.txt

function varargout = inputeqn(varargin)
% INPUTEQN MATLAB code for inputeqn.fig
%      INPUTEQN, by itself, creates a new INPUTEQN or raises the existing
%      singleton*.
%
%      H = INPUTEQN returns the handle to a new INPUTEQN or the handle to
%      the existing singleton*.
%
%      INPUTEQN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INPUTEQN.M with the given input arguments.
%
%      INPUTEQN('Property','Value',...) creates a new INPUTEQN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inputeqn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inputeqn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inputeqn

% Last Modified by GUIDE v2.5 07-Apr-2017 13:01:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inputeqn_OpeningFcn, ...
                   'gui_OutputFcn',  @inputeqn_OutputFcn, ...
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
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.
if isfield(handles, 'data') && ~isreset
    return;
end
handles.data.eqn = '';

% Update handles structure
guidata(handles.figure1, handles);

% --- Executes just before inputeqn is made visible.
function inputeqn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inputeqn (see VARARGIN)

xxdot = evalin('base','xxdot');

% Choose default command line output for inputeqn
handles.output = hObject;

set(handles.text2,'String',['Enter rhs for ' xxdot ' :']);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inputeqn wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inputeqn_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% eqn eingabe
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: 
handles.data.eqn = get(hObject,'String'); % returns contents of edit1 as text
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

%% OK button
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fileID = fopen('eqn.txt','a');              % write data to eqn.txt
fprintf(fileID,'%s\n',handles.data.eqn);    % "
fclose(fileID);
close all;
