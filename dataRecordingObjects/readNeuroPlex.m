%% readNeuroPlex reads .da data files created by the NeuroPlex software from
%  RedShirtImaging. The function is based on the information provided on the
%  RedShirtImaging webpage (http://www.redshirtimaging.com/support/dfo.html)
%  Note: x and y dimensions are exchanged in this file, to ensure identical
%  representations in NeuroPlex and Matlab
%  Syntax:
%  [imageData, bncData, par] = readNeuroPlex (filename)
%
%  If no filename is provided, the function calls uigetfile.
%  par is a structure containing the following information:
%  par.dim       - Size of data set [x y t]
%  par.framedt   - Frame interval
%  par.acqRatio  - Aquisition ratio (electrical over optical)
%  par.darkFrame - Dark frame automatically acquired before the actual
%                  acquistion with closed shutter
% _____________________________________________________________________
% 07/01/2010, Stephan Junek (sjunek@gwdg.de)
function [imageData, bncData, par] = readNeuroPlex (filename)

saveDataToFileMode=1;
%% Select and open file
if nargin ==0
  [filename, pathname] = uigetfile( '*.da', ...
  'Choose a file');
  if isequal(filename,0);return;end
  filename = [pathname filename];  
end

fid = fopen ( filename, 'r', 'l');

if fid == -1
    return;
end

%% Read header info

% Number of columns
fseek (fid, 2*384, 'bof');
par.dim(1) = double(fread(fid,1,'*int16'));

% Number of rows
fseek (fid, 2*385, 'bof');
par.dim(2) = double(fread(fid,1,'*int16'));

% Number of frames (time points)
fseek (fid, 2*4, 'bof');
par.dim(3) = double(fread(fid,1,'*int16'));

% Frame interval (ms)
fseek (fid, 2*388, 'bof');
par.framedt = double(fread(fid,1,'*int16')) / 1000;
fseek (fid, 2*390, 'bof');
dividingFactor = double(fread(fid,1,'*int16'));
if par.framedt >= 10
    par.framedt = par.framedt * dividingFactor;
end

% Aquisition ratio (electrical over optical)
fseek (fid, 2*391, 'bof');
par.acqRatio = double(fread(fid,1,'*int16'));
if par.acqRatio == 0
    par.acqRatio = 1;
end
%% Read image data
fseek (fid, 5120, 'bof');    % goto start offset of par structure
imageData = fread(fid, prod(par.dim),'*int16');

% Data are stored with the dimension order (t,x,y)
if saveDataToFileMode
    imageData = reshape(imageData,par.dim);
    imageData=imageData([1:2:255 256:-2:2],:,:);
else
    imageData = reshape(reshape(imageData, par.dim(3), [])',par.dim);
end


% Exchange x and y to have identical representation as in the NeuroPlex
% software
imageData = permute(imageData, [2 1 3]);
par.dim   = par.dim([2 1 3]);
%% Read BNC data (8 channels)
nrBNC = 8;
bncData = reshape(fread(fid, par.acqRatio*nrBNC*par.dim(3),'*int16'),[par.acqRatio*par.dim(3) nrBNC]);

%% Read dark frame
par.darkFrame = fread(fid, par.dim(1)*par.dim(2),'*int16');
par.darkFrame = reshape(par.darkFrame, par.dim(1:2));

fclose(fid);