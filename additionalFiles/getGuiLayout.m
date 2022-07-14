function success = getGuiLayout( outFolder )
%GETGUILAYOUT Download the 'Nifti and Analyze tools' package
%   This function downloads and add to path the functions inside the
%   fileexchange package 'GUI Layout Toolbox'. It extracts
%   the files from the downloaded zip file and saves them in the folder
%   provided in the arguments. If no folder is specified, it uses the
%   current folder. The output is a logical value set to true if all the
%   steps (download, unzip and save) are successful.
%   For more info about the package, check the page on Mathworks.com:
%   http://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox
%   INPUT: outFolder the folder where the .m files will be stored. If the
%   folder does not exit it will be created (and a warning will be issued)
%   If no input is provided, the package will be downloaded in the current
%   folder.
%   OUTPUT: success a boolean set to true if the download was successful.

%% PRELIMINARY STUFF
success = false;
if nargin < 1 
  saveIn = '.';
elseif ~ischar(outFolder)
  disp('ERROR: outFolder must be a string.')
  return
else
  saveIn = outFolder;
  if ~exist(saveIn, 'dir')
    mkdir(saveIn)
  end
end

if verLessThan('matlab','8.4')
  warning('The most recent version of the GUI only works with R2014b or newer releases')
  disp('Trying to download older version... cannot guarantee functionality!!!')
  remoteZipFile = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/27758/versions/28/download/zip';
else
  remoteZipFile = 'http://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/47982/versions/7/download/zip';
end

%% DOWNLOAD & UNZIP
try
  unzip(remoteZipFile, saveIn)
  % if not in matlab path, add to path
  addpath(genpath(saveIn))
  success = true;
catch
  error('ERROR: An unexpected error occurred while getting the GUI layout.')
end

