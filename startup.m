
% Add the paths to your directories
global baseDirectory
baseDirectory = cd;
addpath(fullfile(baseDirectory, 'Functions'))
addpath(fullfile(baseDirectory, 'Symbolic'))
load_system(fullfile(baseDirectory, 'G2ELib_V1.slx'))

function blkStruct = slblocks
% This function specifies that the library 'mylib'
% should appear in the Library Browser with the 
% name 'My Library'

    Browser.Library = 'G2ELib_V1';
    % 'mylib' is the name of the library

    Browser.Name = 'G2ELib_V1';
    % 'My Library' is the library name that appears
    % in the Library Browser

    blkStruct.Browser = Browser;
end