%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FILE NAME: Video_Player
%   DESCRIPTION: Video Player
%
%   Name          Date          Reason
%   Sahariyaz     28-May-2015   Basic GUI concepts 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = Video_Player(varargin)
% VIDEO_PLAYER MATLAB code for Video_Player.fig
%      VIDEO_PLAYER, by itself, creates a new VIDEO_PLAYER or raises the existing
%      singleton*.
%
%      H = VIDEO_PLAYER returns the handle to a new VIDEO_PLAYER or the handle to
%      the existing singleton*.
%
%      VIDEO_PLAYER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIDEO_PLAYER.M with the given input arguments.
%
%      VIDEO_PLAYER('Property','Value',...) creates a new VIDEO_PLAYER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Video_Player_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Video_Player_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Video_Player

% Last Modified by GUIDE v2.5 28-May-2015 15:38:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Video_Player_OpeningFcn, ...
                   'gui_OutputFcn',  @Video_Player_OutputFcn, ...
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


% --- Executes just before Video_Player is made visible.
function Video_Player_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Video_Player (see VARARGIN)

% Choose default command line output for Video_Player
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Video_Player wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Video_Player_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% ************************ BROWSE *****************************************
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ video_file_name,video_file_path ] = uigetfile({'*.avi'},'Pick a video file');      %;*.png;*.yuv;*.bmp;*.tif'},'Pick a file');
if(video_file_path == 0)
    return;
end
input_video_file = [video_file_path,video_file_name];
set(handles.edit1,'String',input_video_file);
% Acquiring video
videoObject = VideoReader(input_video_file);
% Display first frame
frame_1 = read(videoObject,1);
axes(handles.axes1);
imshow(frame_1);
drawnow;
axis(handles.axes1,'off');
% Display Frame Number
set(handles.text3,'String','1');
set(handles.text4,'String',['  /  ',num2str(videoObject.NumberOfFrames)]);
set(handles.text2,'Visible','on');
set(handles.text3,'Visible','on');
set(handles.text4,'Visible','on');
set(handles.pushbutton2,'Enable','on');
set(handles.pushbutton1,'Enable','off');
%Update handles
handles.videoObject = videoObject;
guidata(hObject,handles);


% ************************ BROWSE TEXT BOX ********************************
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


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

% ************************ START ******************************************
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

videoObject = handles.videoObject;
set(handles.pushbutton3,'Enable','on');
set(handles.pushbutton2,'Enable','off');
axes(handles.axes1);

for frameCount = 2:videoObject.NumberOfFrames
    % Display frames
    set(handles.text3,'String',num2str(frameCount));
    frame = read(videoObject,frameCount);
    imshow(frame);
    drawnow;
end
set(handles.pushbutton4,'Enable','on');
set(handles.pushbutton3,'Enable','off');



% ************************ PAUSE / PLAY ***********************************
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(strcmp(get(handles.pushbutton3,'String'),'Pause'))
    set(handles.pushbutton3,'String','Play');
    uiwait();
else
    set(handles.pushbutton3,'String','Pause');
    uiresume();
end

% ************************ EXIT *******************************************
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);
