function success = checkPackageExist( pkgfile, functionname )
%CHECKPACKAGEEXIST Check if an external package is available and, if not,
%download it
% The function will check for the existance of a function or a class
% provided by a required package. If the file/class is not available, then
% the appropriate doanload function will be called.
% INPUT:
%   pkgfile: a string representing a function, class or file of the
%   required package
%   functionname: the function to call in case the package is not
%   available.
% OUTPUT:
%   success: a boolean value, true if package is available or if it was
%   successfully downloaded, false otherwise.
% USAGE:
%   download.checkPackageExist( 'make_nii', getNiftiTool );
%   download.checkPackageExist( 'specmx.Stack', getSpecMx );
% SEE ALSO:
%   download.getNiftiTool, download.getGuiLayout, download.getSpecMx

success = false;

if nargin ~= 2
  warning('checkPackageExist: 2 input arguments required')
  return;
end

if ~exist(pkgfile, 'class') && ~exist(pkgfile, 'file')
  disp(['You don''t have the toolbox providing ', pkgfile])
  reply = input('Do you want to download it now? Y/N [Y]:','s');
 if isempty(reply)
    reply = 'N';
 end
 if ~strcmpi(reply, 'Y')
   disp('WARNING: Cannot proceed without the package.')
   return;
 else
   outFolder = uigetdir('.', 'Select on which folder the toolbox will be saved');
   if 0 == outFolder
     warning('WARNING: No directory specified, downloading package in current directory.')
     outFolder = '.';
   end
   funcToCall = str2func(strcat('download.', functionname));
   success = funcToCall(outFolder);
   if ~success
     disp('WARNING: Cannot download tool. Aborting now.')
     return
   end
 end
else %nothing to do
  success = true;
end

end

