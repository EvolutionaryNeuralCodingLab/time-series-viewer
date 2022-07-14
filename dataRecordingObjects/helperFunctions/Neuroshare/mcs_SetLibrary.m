function [ns_RESULT] = mcs_SetLibrary(filename);
% Opens a Neuroshare Shared Library (.DLL or .so) in prepearation for other work
%
%   Usage:
%                [ns_RESULT] = mcs_SetLibrary('filename.dll')
%
%   Description:
%                Opens the dynamic linked library specified by filename
%
%   Parameters:
%                filename      Pointer to a null-terminated string that specifies the
%                              name of the file to open. 
%
%   Return Values:
%                ns_RESULT   This function returns 0 if the file is successfully
%                            opened. Otherwise there is an error message shown in 
%                            the CommandWindow and the function tries to open the 
%                            file using different paths.


[ns_RESULT] = ns_SetLibrary(filename);                                     % first step: Opens a Neuroshare shared Library, load the appropriate DLL

[pathpart,namepart,extpart]= fileparts(filename);                          % returns the path name, file name, and extension for the specified file

if (ns_RESULT ~= 0)                                                        % if nsresult is unequal 0, Dll was not found

    disp(['ns_SetLibrary does not find the library: ' filename]);          % error message is shown if nsresult is unequal 0    
    if (strcmp(pathpart, ''))                                              % compares pathpart for equality with an empty string 

        disp('Searching the Matlab path...');                              % as long as pathpart is empty (not yet found),information, that search is in progress
        
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        a = path;                                                          % parse the Matlab path,path displays the MATLAB search path
                                                                           % and returns the search path to string variable a
        jj=1;                                                              % definition of jj
        kk=1;                                                              % definition of kk
        c='';                                                              % definition of c as an empty string
        
        for ii = (1:length(a))                                             % definition of ii (1: length of path)

            if (a(ii) == pathsep)                                          % if condition proofs if 'a'(=path) is a match to pathsep,
                                                                           % pathsep returns the search path separator character for this platform  
                                                                           % The search path separator is the character that separates path names 
                                                                           % in the pathdef.m file, as returned by the path function.
                b{jj} = c;                                                 % definition of b as an empty string
                c = '';                                                    % definition of c as an empty string
                jj = jj + 1;                                               % jj increases by 1
                kk = 1;                                                    % kk stays the same
            
            else                                                           % otherwise
                c(kk) = a(ii);                                             % c(kk) is defined as a(ii)(=path)
                kk = kk + 1;                                               % kk increases by 1
            end
            ii = ii + 1;                                                   % ii increases by 1

        end
%----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        for ii = (1:length(b))                                             % try the different paths
           c = [b{ii} filesep filename];                                   % new vector c containing b, filesep, filename
                                                                           % filesep returns the platform-specific file separator character
                                                                           % The file separator is the character that separates individual folder and file names in a path string.
           
           if (exist(c) == 3)                                              % if file exists try ns_SetLibrary again
                                                                           % output 3 means 'c exists as a MEX- or DLL-file on your MATLAB search path'
               
               
               
               disp(['Try to open ' c]);                                   % information that another try is in progress 
               
               [ns_RESULT] = ns_SetLibrary(c);                             % return to ns_SetLibrary, called by c 
               
               if (ns_RESULT == 0)                                         % if condition proofs if ns_RESULT is 0
                   disp('Success');                                        % 'Success' gives you certainty that library is found  
                   break;
               end
           end
        end
    else
        disp('Using specific path, further searching is skipped');
    end
    
end

end

