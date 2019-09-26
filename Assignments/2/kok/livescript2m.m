function livescript2m(file)
    %LIVESCRIPT2M Convert Live Script file to a classic m-file script
    %
    % LIVESCRIPT2M(FILE) extracts, in the same directory, a MATLAB script FILE.m
    % from the Live Script file FILE.mlx
    %
    % Peter Corke (c) 2016
    
    % an .mlx file is simply a zip archive with the MATLAB code inside an XML
    % file matlab/document.xml.  Use MATLAB's built in XSL engine and an XSLT
    % file to parse out the MATLAB code.
    
    mFile = [file '.m'];
    mlxFile = [file '.mlx'];
    
    if exist(mFile, 'file') == 2
        error('attempting to overwrite existing m-file %s', mFile);
    end
    
    % find the XSLT file, in same dir as this function
    pathstr = fileparts( mfilename('fullpath') );
    xsltFile = fullfile(pathstr, 'livescript2m.xslt');
    
    % create a place to safely unzip the live edit file
    unzipDir = fullfile(tempname);
    
    % unzip it
    unzip(mlxFile, unzipDir)
    
    % and point to the xml code inside the directory structure
    srcFile = fullfile(unzipDir, 'matlab', 'document.xml');
    
    % now we'll transform the xml to matlab
    xslt(srcFile, xsltFile, mFile);
    
    rmdir(unzipDir, 's')
end