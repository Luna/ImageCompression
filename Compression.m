function varargout = Compression(varargin)
% COMPRESSION MATLAB code for Compression.fig
%      COMPRESSION, by itself, creates a new COMPRESSION or raises the existing
%      singleton*.
%
%      H = COMPRESSION returns the handle to a new COMPRESSION or the handle to
%      the existing singleton*.
%
%      COMPRESSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPRESSION.M with the given input arguments.
%
%      COMPRESSION('Property','Value',...) creates a new COMPRESSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Compression_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Compression_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Compression

% Last Modified by GUIDE v2.5 05-Dec-2022 23:38:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Compression_OpeningFcn, ...
    'gui_OutputFcn',  @Compression_OutputFcn, ...
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


% --- Executes just before Compression is made visible.
function Compression_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Compression (see VARARGIN)

% Choose default command line output for Compression
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Compression wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Compression_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Upload Images
% --- Executes on button press in uploadButton.
function uploadButton_Callback(hObject, eventdata, handles)
global UncompressedImage

% opens a file and then displays the original image.
[file,path] = uigetfile({'*.bmp';'*.png';'*.jpeg';'*.*'},...
    'Select an Image file');

% User selected Cancel, return to the previous screen and display error
if isequal(file,0)
    f = errordlg('File not found','File Error');
else
    filepath = fullfile(path,file);
    disp(['User selected ', fullfile(path,file)]);
end

OriginalImage = imread(filepath);
UncompressedImage = rgb2gray(OriginalImage); %convert this RGB image into grayscale
axes(handles.originalImage); % Displays Gray Image on Figure Axis
imshow(UncompressedImage); title('Original Image');
imwrite(UncompressedImage,"GrayStart.png"); % Writes gray image to disk, used as original image from this point

%% Fourier Compression
% --- Executes on button press in fourierCompression.

function fourierCompression_Callback(hObject, eventdata, handles)
global UncompressedImage;
% FFT analysis on the uncompressed image
At=fft2(UncompressedImage);

% The Fourier coefficients are then zeroed out and inverse transformed to make them small enough for compression.
thresh=.1*0.005* max(max(abs(At)));
ind = abs(At)>thresh; % values greater than 0.1*0.005 are kept in the compressed image
Atlow = At.*ind; % Calculates number of small coefficients
Alow=uint8(ifft2(Atlow)); % compresses the original image into a smaller one using an IFFT algorithm
axes(handles.FFT); %
imshow(Alow); % Displays compressed image to FFT Axis
title('FFT Compression');
imwrite(Alow,"fftCompression.png"); % Writes compressed Image to the disk

% Used for calculation of metrics
z = dir("fftCompression.png");
fftSize=z.bytes;


US=dir("GrayStart.png");
UncompressedSize = US.bytes;


% Compression Ratio - ratio of the size of the compressed image divided by the size of the uncompressed image
CR =  UncompressedSize / fftSize ;
h = msgbox(sprintf('Fourier Compression Ratio = %.5f', CR));

% MSE
P = im2double(UncompressedImage);
Alow = im2double(Alow);
err = immse(Alow, P);
h = msgbox(sprintf('Fourier Compression MSE = %.4f', err));

% PSNR
[peaksnr, snr] = psnr(Alow, P);
h = msgbox(sprintf('Fourier Compression PSNR = %.4f', peaksnr));

%FourierMetrics = [CR, err, peaksnr, snr];


%% SVD Compression
% --- Executes on button press in svdCompression.

function svdCompression_Callback(hObject, eventdata, handles)
global UncompressedImage
doubleUncompressedImage = im2double(UncompressedImage); % uncompressed image and returns a double-precision floating point value
[u, s, v] = svd(doubleUncompressedImage); %  The svd function then computes the singular values of this image
sz = size(s) * 0.25; % Scaled down to 25% of its original size (sz)

s2=s;
s2(sz:end, :) = 0;
s2(:, sz:end) = 0;
SVDCompressed=u*s2*v';
axes(handles.SVD);
imshow(SVDCompressed);
title('SVD Compression');
imwrite(SVDCompressed,"SVDCompression.png")

% Used for calculation of metrics
z = dir("SVDCompression.png");
fftSize=z.bytes;


US=dir("GrayStart.png");
UncompressedSize = US.bytes;


% Compression Ratio
CR =  UncompressedSize / fftSize;
h = msgbox(sprintf('SVD Compression Ratio = %.5f', CR));

% MSE
P = im2double(UncompressedImage);
SVDCompressed = im2double(SVDCompressed);
err = immse(SVDCompressed, P);
h = msgbox(sprintf('SVD Compression MSE = %.4f', err));

% SNR
[peaksnr, snr] = psnr(SVDCompressed, P);
h = msgbox(sprintf('SVD Compression SNR = %.4f', snr));

% PSNR
h = msgbox(sprintf('SVD Compression PSNR = %.4f', peaksnr));

%% RLE Compression
% --- Executes on button press in rleCompression.

function rleCompression_Callback(hObject, eventdata, handles)
global UncompressedImage
str=reshape(UncompressedImage.',1,[] );
str=reshape(UncompressedImage.',1,[] ); % Convert Image from 2D to 1D

y=[];
c=1;
for i=1:length(str)-1
    if(str(i)==str(i+1))
        c=c+1;
    else
        y=[y,c,str(i),];
        c=1;
    end
end
y=[y,c,str(length(str))];
writematrix(y,'RLE.txt')
f = msgbox('Run Length Encoding Done');
ff = size(str,[2]);
gg = size(y, [2]);
cr = ff/gg

h = msgbox(sprintf('RLE Compression Ratio = %.5f', cr));
[fx,fy] = size(UncompressedImage);
save('RLE_RES.mat','fx', 'fy');

%% RLE Decompression
% --- Executes on button press in decompressButton.

function decompressButton_Callback(hObject, eventdata, handles)
[file,path] = uigetfile('*.txt',...
    'Select Compressed file');
if isequal(file,0)
    f = errordlg('File not found','File Error');
else
    filepath = fullfile(path,file);
    disp(['User selected ', fullfile(path,file)]);
end

y = readmatrix(filepath);
b=[];
for i=1:2:length(y)
    b=[b y(i)];
end
z=[];
for i=2:2:length(y)
    z=[z y(i)];
end
v=[];
parfor i=1:length(b)
    m=b(i);
    q=z(i);
    for j=1:m
        v=[v q];
    end
end
str=v;
DD=uint8(str);
load('RLE_RES.mat','fx')
load('RLE_RES.mat','fy')
Restore=reshape(DD,fy,fx);
axes(handles.largeOriginal);
axis off;
imshow(Restore.');
title('Decoded Image');
imwrite(Restore,"Decompressed.png")
f = msgbox('Image has been decompressed');

% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
clear all; close all; clc
Compression
