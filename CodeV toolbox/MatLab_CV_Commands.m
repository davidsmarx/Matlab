%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: MatLab_CV_Commands.m 
%
% Purpose:
% Demonstrate usage of the MatLab - CODE V interface.  The program
% illustrates many commands that can be used to exchange data between
% MatLab and CODE V and how MatLab can issue commands to CODE V.
% 
% Commands Demonstrated:
%
% Utility Functions
%   Start/StopCodeV
%   Get/SetCommandTimeout
%   Get/SetMaxTextBufferSize
%   Get/SetStartingDirectory
%   GetCodeVVersion
%
% Synchronous Commands
%   Command
%   EvaluateExpression
%   GetCurrentOption
%   GetCurrentSubOption
%   GetZoomCount
%   GetSurfaceCount
%   GetFieldCount
%   GetWavelengthCount
%   GetDimension
%   GetStopSurface
%   GetMaxAperture
%   BufferToArray
%   ArrayToBuffer
%
% Asynchronous Commands
%   AsynCommand
%   IsExecuting
%   Wait
%   GetCommandOutput
%   StopCommand
%
% The following Macro Functions
%   ZFRFIT/ZRNFIT 
%   GAUSSWTS
%   BESTSPH
%   EVALZERN
%   FITERROR
%   NORMRADIUS
%   RAYRSI
%   RAYSIN
%   RAYTRA
%   ZERNIKE
%   ZFRCOEF
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables/Programs
%
%	An input parameter is passed to the program.  An output parameter is
% 	passed back or changed by this program for use by another.
%
%	Inputs: None	
%	Outputs: None
%	Called programs: CODE V
%	Global Variables: None 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by: Byron Taylor
% Incept date: 3 Feb 2003
% Status: Complete
% Modifications:
%   26 Feb. 2003: COM interface was changed to maintain consistant
%   Row/Column addressing as opposed to the previous case in which rows and
%   columns were transposed
%
% Notes:
%   1) Not all possible commands have been included.  See CV User
%   Command.doc for the full list.
%   2) Will need to enter a path for CODE V working directory, i.e. lens
%   file locations for some examples -- see "path" variable below
%   3) Will need to load a lens/sequence file for some examples -- see
%   "lens_to_load" variable below.
%   4) To activate commands, find the command in the body and select it
%   by setting the "if" argument to "1"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create ActiveX Server
%What is CVCommand
%CVCommand is a COM interface on CODE V which is designed to allow command mode access 
%to CODE V functionality.  It can be used to automate tasks or to get data which is 
%used in calculations in other programs.  Each time a session of CVCommand is created, 
%it will start a new copy of CODE V on your machine (this is done when the StartCodeV 
%function is called).  CVCommand has no GUI and therefore no graphics support any plots 
%which are created must be output to a file or they are lost.  Plot files can be viewed 
%either in a session of the CODE V GUI or in the CODE V plot viewer.

%When you call actxserver, you request an instance of the CVCommand interface.  This 
%instance should be requested by object name (e.g. "CodeV.Command.9XX").  The version 
%number corresponds to the CODE V version.  A specific version number should be 
%requested because the interfaces may change with later versions of CODE V. 

%Once an instance is created, the utility functions can be called on it to setup the 
%environment (timeout, buffer size, starting directory, etc.).  Then, CODE V can be 
%started using the StartCodeV function.  It can be stopped using the StopCodeV 
%function, which should be done before your program exits to ensure a clean shutdown 
%of the CODE V session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paths = 'c:\cvuser' ;                     % will need to load a path for files
lens_to_load = 'cv_lens:dbgauss' ;        % will need to load a lens for some examples

cv=actxserver('CODEV.Command') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Get/SetCommandTimeout
%
% DESCRIPTION
% These functions are used to get or set the timeout for synchronous commands.  
% They have no effect on asynchronous commands.  Can be set before CODE V
% is started.
%
% timeout_set       Input: Set timeout in seconds ??? (integer)
% Result            Output: What the timeout is set to (integer)  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)

    timeout_set = 500000 ;                              

    invoke(cv,'SetCommandTimeout',timeout_set) ;        
    [Result] = invoke(cv,'GetCommandTimeout') ; 

     disp(' ')
     disp(['Set Timeout to = ' int2str(timeout_set)])
     disp(['Read Timeout Set to = ' int2str(timeout_set)])

    clear timeout_set Result
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get/SetMaxTextBufferSize
%
% DESCRIPTION
% These functions are used to get or set the maximum buffer size for text returned 
% by the Command and GetCommandOutput functions. Can be set before CODE V
% is started.
%
% buffersize_set      Input: Buffer size in characters, default is 256000. (Long Integer)
% Result              Output: What buffer size is set to (Long Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    buffersize_set = 300000 ;     

    invoke(cv,'SetMaxTextBufferSize',buffersize_set) ;            
    [Result] = invoke(cv,'GetMaxTextBufferSize') ; 

     disp(' ')
     disp(['Set Buffer Size to = ' int2str(buffersize_set)])
     disp(['Read Buffer Size Set to = ' int2str(Result)])
     
%    disp(['Set Buffer Size to = ' int2str(buffersize_set)])
%    disp(buffersize_set)
%    disp(['Read Buffer Size Set to = ' Result])
%    disp(Result) 
  
    clear buffersize_set Result
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Get/SetStartingDirectory
%
% DESCRIPTION
% Get or set the working directory for CODE V.  This is most often called before 
% StartCodeV to set the directory of execution. Can be set before CODE V
% is started.
%
% set_starting_directory    Input: Directory -- Input (String)
% Result                    Output: Returned Starting Directory (String)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    set_starting_directory = paths ;             
    
    clear paths ;
    
    invoke(cv,'SetStartingDirectory',set_starting_directory) ;  % Input is String
    
    [Result]= invoke(cv,'GetStartingDirectory') ;  % Input/Output are strings

    disp(' ')
    disp(['Starting Directory Set to = ' Result])

    clear set_starting_directory Result ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Start CODE V
%
% DESCRIPTION
% These functions start or stop the CODE V session being run by CVCommand.  
% Start must be called before any function other than Set/GetCommandTimeout, 
% Set/GetMaxTextBufferSize, GetCodeVVersion, or Set/GetStartingDirectory is called.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

invoke(cv,'StartCodeV') ;
invoke(cv,'command',['res ' lens_to_load])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GetCodeVVersion
%
% DESCRIPTION
% This function is used to get the version of CVCommand.  This version number is 
% identical to the version of CODE V that it is running.
%
% Result    Output: Current CODE V version (string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    [Result] = invoke(cv,'GetCodeVVersion') ; 
    
    disp(' ')
    disp(['CODE V Version = ' Result])

    clear Result ;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synchronous Usage Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Command
%
% DESCRIPTION
% This function sends a command to the CODE V session being run by CVCommand and 
% returns its output.  Calling this function clears the results of the previous 
% asynchronous function call.
%
% Input     Input: Any valid CODE V command, except plotting. (String)
% Result    Output: Echo back to CODE V command window. (String)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    Input = 'sur sa ' ;
    [Result] = invoke(cv,'command',Input) ;

    disp(' ')
    disp(Input)
    disp(Result)

    clear Result Input ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EvaluateExpression
%
% DESCRIPTION
%
% This evaluates an expression and returns its value.  It is equivalent to the 
% EVA command in CODE V. Note that because this is a string, the value is only as 
% precise as the output into the string - it is not a true floating point value.
%
% Input     Input: Any valid CODE V "eva" command. (String)
% Result    Output: Echo back to CODE V command window. (String)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    Input = '(op s1..i-2)' ;
    [Result] = invoke(cv,'EvaluateExpression',Input)  ; 

    disp(' ')
    disp([Input ' = ' Result])
    
    clear Result Input ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CODE V State Information
% 
% All of these functions will fail if CODE V is currently running an asynchronous 
% command.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% GetCurrentOption
%
% DESCRIPTION
% Returns the name of the current option. Returns CHA in not currently in an
% option
%
% Result    Output: Cuurent Option. CHA if not in one (String)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1) 
    
    invoke(cv,'command','PMA') 

    clear Result ;
    
    [Result] = invoke(cv,'GetCurrentOption') %Output is string
    invoke(cv,'command','can') ;
    
    disp(' ') 
    disp(['Current option selection = ' Result])

    clear Result ;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GetCurrentSubOption
%
% DESCRIPTION
% Returns the name of the current suboption. Returns an empty string if 
% CODE V is not currently in a sub-option.
%
% Result    Output: Cuurent Option. Empty if not in one (String)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [Result] = invoke(cv,'GetCurrentSubOption') ; % Output is string

    disp(' ') 
    disp(['Current suboption selection = ' Result])
    
    clear Result ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% System Parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    clear parameters ;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% GetFieldCount
%
% DESCRIPTION
% Returns the number of fields
%
% fields:   Output: Returns the current number of fields. (Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fields = invoke(cv,'GetFieldCount') ; 
    parameters(1,1) = {'Number Fields'} ;
    parameters(1,2) = {fields} ;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% GetWavelengthCount
%
% DESCRIPTION
% Returns the current number of wavelengths.
%
% wave:     Output: Returns the current number of wavelength. (Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    waves = invoke(cv,'GetWavelengthCount') ; 
    parameters(2,1) = {'Number Wavelengths'} ; 
    parameters(2,2) = {waves} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GetDimension
% 
% DESCRIPTION
% Returns a value representing the type of dimensions in the system (integer)
%
% Dimension:    Output (integer)
%                   0.	inches
%                   1.	centimeters
%                   2.	millimeters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dimension = invoke(cv,'GetDimension') ; 
    switch dimension
        case 0
            dims = 'Inches' ;
        case 1
            dims = 'Centimeters' ;
        case 2
            dims = 'Millimeters'; 
        otherwise
            error('Unexpected units')
    end
    parameters(3,1) = {'Dimensions'} ;
    parameters(3,2) = {dims}   ;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GetStopSurface
% 
% DESCRIPTION
% Returns the surface number of the current stop surface.
%
% stop_sur:     Output: Returns the stop surface number. (integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    stop_sur = invoke(cv,'GetStopSurface') ; 
    parameters(4,1) = {'Stop Surface Number'} ;
    parameters(4,2) = {stop_sur} ;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GetZoomCount
%
% DESCRIPTION
% Returns the current number of zoom positions.
%
% zooms:     Output: Returns the number of zooms (integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    zooms = invoke(cv,'GetZoomCount') ; 
    parameters(5,1) = {'Number Zooms'} ;
    parameters(5,2) = {zooms}  ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GetSurfaceCount
%
% DESCRIPTION
% Returns the current number of surfaces. 
%
% num_surfaces:     Output: Number of surfaces (same as (num s)). (Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    num_surfaces = invoke(cv,'GetSurfaceCount') ;
    parameters(6,1) = {'Number Surfaces'} ;
    parameters(6,2) = {num_surfaces} ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GetMaxAperture
%
% DESCRIPTION
% Returns the maximum aperture size for the specified surface and zoom.
% Surface	The number of the surface that the maximum aperture will be determined for.
% Zoom	The Zoom position at which the maximum aperture will be determined.
% MaxAperture	Contains the maximum aperture size.  This uses the MAP database item.
% returns aperture radius in lens units.
%
% Surface:      Input: Surface number to compute MAP on (Integer)
% Zoom:         Input: Zoom number to use (Integer)
% Result:       Output: Radius of aperture (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Zoom = zooms ;
    invoke(cv,'command','lis') ;        % Need this command for some reason before
                                         % maximum aperture will run correctly
    for Surface=1:num_surfaces-1
        MAP{Surface,1} = ['Surface S' int2str(Surface)] ;
        MAP{Surface,2} = invoke(cv,'GetMaxAperture',Surface,Zoom) ; 
    end

    disp(' ')
    disp('Maximum Aperture Radii')
    disp(MAP) 
    disp('System Parameters')
    disp(parameters)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    clear MaxAperture_size Zoom Surface MAP ii num_surfaces stop_sur dimension ...
                waves fields zooms parameters dims;

end               
%%%%%%%%%%%%% End of System Parameters Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ArrayToBuffer
%
% DESCRIPTION
% This function copies the numeric contents of an array to a CODE V worksheet buffer.
% Rows	    Input: Number of rows to copy to the buffer (Integer)
% Cols	    Input: Number of colums to copy to the buffer (Integer)
% Buffer	Input: Worksheet buffer number to be copied (Integer)
% Array	    Input: Array of doubles to be fill the worksheet.  This array must 
%           be (nEndRow - nStartRow + 1)x(nEndCol - nStartCol + 1).
% Return    Output: Contains a value indicating success or failure: (Integer)
%               0	success
%               -1	failure: the buffer had zero rows or columns
%               n	The buffer was not empty and the array was concatenated with 
%                   the buffer after row n of the buffer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    Buffer = 117 ;   % Buffer to put data in
    
    SelectArray = 'General' ; %Square (NxN), Vector1 (1xN), Vector2 (Nx1), General (MxN)
    
    switch SelectArray
        
        case 'Square' % Pass a NxN array to CODE V
            
            Dim = 4 ;
            Rows = Dim ;
            Cols = Dim ;
            kk = 0 ;
            for ii = 1:Dim
                for jj = 1:Dim
                    Array(ii,jj) = kk ;
                    kk = kk + 1 ;
                end
            end
             
        case 'Vector1' % Pass a 1xN vector to CODE V 
            Dim = 4 ;
            Rows = 1 ;
            Cols = Dim ;
            Array = [1:1:Dim] ;
            
        case 'Vector2' % Pass a Nx1 vector to CODE V 
            Dim = 4 ;
            Rows = Dim ;
            Cols = 1 ;
            Array = [1:1:Dim]' ;
            
        case 'General' % Pass an array MxN to CODE V
            Rows = 4 ;
            Cols = 2 ;
            kk = 0 ;
            for ii = 1:Rows
                for jj = 1:Cols
                    Array(ii,jj) = kk ;
                    kk = kk + 1 ;
                end
            end
   
        otherwise
            error('Invalid SelectArray choice')
    end
 
   invoke(cv,'command',['buf del b' int2str(Buffer)]) ;
   
   clear result ;
   
   [result] = invoke(cv,'ArrayToBuffer',Rows ,Cols ,Buffer,Array);

   disp(' ')
   disp('Return =')
   disp(result)
   disp('Passed Array')
   disp(Array) 
   disp('Buffer Result')
   invoke(cv,'command',['buf lis b' int2str(Buffer)]) 

   clear Buffer Dim Rows Cols Buffer Array kk jj ii SelectArray ;
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BufferToArray
%
% DESCRIPTION
% This function copies the numeric contents of a CODE V worksheet buffer to an 
% array.  Buffer Rows map to the array's second subscript. Buffer Columns map to 
% the array's first subscript.
%
% StartRow	Input: Starting row number of the buffer (Integer)
% EndRow	Input: Ending row number of the buffer (Integer)
% StartCol	Input: Starting column number of the buffer (Integer)
% EndCol	Input: Ending column number of the buffer (Integer)
% Buffer	Input: Worksheet buffer number to be copied (Integer)
% Array	    Output: Array of doubles to be filled with the worksheet data.  
%           This array must be (EndRow - StartRow + 1)x(EndCol - StartCol + 1).
% Return 	Output: 0 for success, 1 for failure (usually because of
%           non-numeric data in the buffer).  (Integer)
%
% In this example, it would be more compact to use the ArrayToBuffer
% command to fill the buffer, but this illustrates some command usage.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    Buffer = 117 ;   % Buffer to get data from
    
    invoke(cv,'command',['buf del b' int2str(Buffer)]) ;
    
    SelectArray = 'General' ; %Square (NxN),Vector1 (1xN),Vector2 (Nx1),General (MxN)
    
    switch SelectArray
        
        case 'Square' % Load a NxN buffer in CODE V
            
            Dim = 4 ;
            Rows = Dim ;
            Cols = Dim ;
            kk = 0 ;
            for ii = 1:Dim
                for jj = 1:Dim
                    invoke(cv,'command',['buf put b' int2str(Buffer) ' i' ...
                                     int2str(ii) ' j' int2str(jj) ' ' int2str(kk)]) ;
                    kk = kk + 1 ;            
                end
            end
             
        case 'Vector1' % Load a 1xN buffer in CODE V
            Dim = 4 ;
            Rows = 1 ;
            Cols = Dim ;
            kk = 0 ;
            for ii = 1:Dim
                invoke(cv,'command',['buf put b' int2str(Buffer) ' i' ...
                                     int2str(1) ' j' int2str(ii) ' ' int2str(kk)]); 
                kk = kk + 1 ;    
             end
            
        case 'Vector2' % Load a Nx1 buffer in CODE V 
            Dim = 4 ;
            Rows = Dim ;
            Cols = 1 ;
            kk = 0 ;
            for ii = 1:Dim
                invoke(cv,'command',['buf put b' int2str(Buffer) ' i' ...
                                    int2str(ii) ' j' int2str(1) ' ' int2str(kk)]); 
                kk = kk + 1 ;
             end
                   
        case 'General' % Pass an array MxN to CODE V
            Rows = 4 ;
            Cols = 2 ;
            kk = 0 ;
           for ii = 1:Rows
                for jj = 1:Cols
                    invoke(cv,'command',['buf put b' int2str(Buffer) ' i' ...
                                        int2str(ii) ' j' int2str(jj) ' ' int2str(kk)]); 
                kk = kk + 1 ;
                end
            end
   
        otherwise
            error('Invalid SelectArray choice')
    end
 
   StrtRow = 1;
   StrtCol = 1;
   EndRow = Rows ;
   EndCol = Cols ;
   
   Array = zeros(Rows,Cols) ;
   
   clear return_val ;
   
   [return_val, Array] = ...
        invoke(cv, 'BufferToArray', StrtRow ,EndRow ,StrtCol ,EndCol ,Buffer, Array) ;
   
   disp(' ') 
   disp('Return = ')
   disp(return_val)    
   disp('Buffer Data')
   invoke(cv,'command',['buf lis b' int2str(Buffer)]) 
   disp('Array Read In')
   disp(Array) 
   
   clear Buffer Dim Rows Cols Buffer Array kk jj ii SelectArray StrtRow StrtCol ...
                                                   return_val EndRow EndCol ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Asynchronous Commands
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% AsyncCommand
% 
% DESCRIPTION
% Start an asynchronous command.  Only one can be run per session of CVCommand; 
% this function call fails if CODE V is already running a command.  Calling this 
% function clears the results of the previous asynchronous function call.
%
% CommandLine   Input: Command to issue to CODE V to be run asynchronously
%               (string)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)                         
    
    invoke(cv,'command','res cv_lens:dbgauss') ; 
    invoke(cv,'command','var sa');
    
    CommandLine = 'aut;mnc 2;go' ;     
    invoke(cv,'AsyncCommand',CommandLine) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Wait
%
% DESCRIPTION
% Wait for an asynchronous command to complete.
%
% WaitTime      Input: The time to wait in seconds ??
% Result        Output: (0) means command completed
%                       (1) means command still running
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear result ;
    WaitTime = 2 ;
    [result] = invoke(cv,'Wait',WaitTime) ;
    disp(['WaitTime result = ' num2str(double(result))])
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% IsExecutingCommand
%
% DESCRIPTION
% This function returns a boolean value indicating whether an asynchronous command 
% is currently executing.
%
% Result        Output: Boolean which indicates whether or not an asynchronous command is 
%               currently executing.
%                   0 Finished
%                   1 Still running
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear result            
    
    [result] = invoke(cv,'IsExecutingCommand')  ; 
    count = 0 ;
    while [result] ~= 0
        [result] = invoke(cv,'IsExecutingCommand');
        count = count + 1 ;
    end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
% GetCommandOutput
%
% DESCRIPTION
% Gets the output from a completed command.  This function is normally used 
% to get the output from an asynchronous command run.  However, note that it will 
% return the output of the last completed command.  Therefore, you should not call 
% any function which might reset this buffer, such as Command, EvaluateExpression, 
% or any of the math and optical functions.
%
% Result    Output: Contains the output.(String) Maximum length is the maximum 
%           buffer size.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    clear result ;
    [result] = invoke(cv,'GetCommandOutput') ;
    result

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% 
% StopCommand
%
% DESCRIPTION
% This function aborts the currently running CODE V calculation. No return
% value.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    invoke(cv,'StopCommand') ;
    
    disp(['Number of times checked (Total run time is WaitTime*Count)) ' ...
                                                        int2str(count)])
    
    clear CommandLine WaitTime result count
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Math and Optical MACRO Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Zernike Fitting Functions: ZRNFIT/ZFRFIT
%
% DESCRIPTION.  
% This is equivalent to calling the ZFRFIT or ZRNFIT macro functions in CODE V.
%
% NumPoints	    Input: The number of points that will be fitted (Integer)
% X	            Input: An vector of x values NumPoints long (Double)
% Y	            Input: An vector of y values NumPoints long (Double)
% F	            Input: An vector of f values NumPoints long (Double)
% NumZTerms     Input: Number of zernike coefficients to fit to.  Must be dimensioned
%               long enough to hold correct number of terms. CODE V will return enough
%               coefficients to complete an order and will round up if needed. (Integer) 
% Coeffs	    Output: An array of coefficients NumZTerms long. (Double)   
% Result	    Output: The result of the CODEV macro function call.  The return value 
%               here is the RMS fit error; if the fit fails, the return value is set 
%               to -1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    NumPoints = 100 ;
    scale = 1/2^0.5 ;
    x = 2*scale*rand([NumPoints 1]) - scale ;
    y = 2*scale*rand([NumPoints 1]) - scale ;

    f = zeros(NumPoints,1) ;

    zernike_type = 'defocus' ;

    switch zernike_type
        case 'piston'
            f(:) = 0.11 ;       % Use to test piston
            NumZRNTerms = 3 ;   % Number of ZRN terms to fit --  1 is not a valid choice
            NumZFRTerms = 4 ;   % Number of ZFR terms to fit --  1 is not a valid choice    
        case 'x tilt'
            r = (x.^2 + y.^2).^0.5 ;
            cos_theta = x./r  ;     %Assumption is that r will never be exactly zero
            f = r.*cos_theta ;
            NumZRNTerms = 3 ;       % Number of ZRN terms to fit
            NumZFRTerms = 4 ;       % NNumber of ZRN terms to fit 
        case 'defocus'    
            r = (x.^2 + y.^2) ;
            f = 2.*r - 1 ;
            NumZRNTerms = 6 ;       % Number of ZRN terms to fit
            NumZFRTerms = 4 ;       % Number of ZRN terms to fit        
        otherwise
        error('Invalid Zernike Term Choice') ;
    end

    Coeffs = zeros(NumZRNTerms,1) ;
    disp(['ZRNFIT results for ' zernike_type])
    clear result ;
    [result,Coeffs] = invoke(cv,'ZRNFIT',NumPoints, x, y, f, NumZRNTerms, Coeffs)  ;

    for ii = 1:NumZRNTerms
        C(ii,1) = {['Z' int2str(ii)]} ;
        C(ii,2) = {Coeffs(ii)} ;
    end

    disp(' ')
    disp('ZRN Coefficients = ')
    disp(C) 
    
    Coeffs = zeros(NumZFRTerms,1) ;
    disp(['ZFRFIT results for ' zernike_type])
    clear result ;
    [result,Coeffs] = invoke(cv,'ZFRFIT',NumPoints, x, y, f, NumZFRTerms, Coeffs) ;

    clear C ;
    for ii = 1:NumZFRTerms
        C(ii,1) = {['Z' int2str(ii)]} ;
        C(ii,2) = {Coeffs(ii)} ;
    end

    disp(' ')
    disp('ZFR Coefficients = ')
    disp(C) 
    
    clear x y f r cos_theta scale NumZRNTerms NumZFRTerms NumPoints ...
                                                      C zernike_type Coeffs result
    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GAUSSWTS
%
% DESCRIPTION
% This is equivalent to calling the GAUSSWTS macro function in CODE V.
%
% NumInputPts	Input: The number of input points (Integer)
% InputCoords	Input: A NumInputPts vector of coordinates at which weights are   
%               supplied.  The input coordinates need not be equally spaced. (Double)
% InputWeights 	Input: An NumInputPts vector of weights at the given input 
%               coordinates. (Double)
% NumQuadPts	Input: The number of Gaussian quadrature points and weights 
%               desired.(Integer)
% Coords	    An output vector, dimensioned at NumQuadPts.  It will receive the 
%               coordinates to be used for the numerical integration. (Double)
% Weights 	    An output array, dimensioned at NumQuadPts.  It will receive the 
%               weights to be used for the numerical integration. (Double)
% Result	    The return value of the GAUSSWTS macro function.  It is 0 if there 
%               are no erros in the computation and -1 if any errors are 
%               encountered. (Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)  
    
    NumInputPts = 11 ;
    NumQuadPts = 11 ;
    
    lower_bound = -1 ;
    upper_bound = 1 ;
    step = (upper_bound-lower_bound)/(NumInputPts - 1) ;
    %Inputs
    InputCoords = (lower_bound:step:upper_bound)' ;
    InputWts = ones(NumInputPts, 1) ;  % Could be any function, but w(x) = 1 was chosen
    %Outputs
    Coords = zeros(NumQuadPts,1) ;
    Weights = zeros(NumQuadPts,1) ;

    clear Result ;
    [Result,Coords,Weights] = invoke(cv,'GAUSSWTS',NumInputPts,InputCoords, ...
                                                 InputWts,NumQuadPts,Coords,Weights) ; 

    C(1,1) = {'Coordinate'} ;
    C(1,2) = {'Weights'} ;
    for ii = 1:NumQuadPts
        C(ii+1,1) = {(Coords(ii))} ;
        C(ii+1,2) = {(Weights(ii))} ;
    end
      
    disp('Gaussian Weights')
    disp(C)
    
    clear NumInputPts NumQuadPts lower_bound upper_bound step InputCoords ...
             ii InputWts Coords Weights Result Coords Weights C ;
end
                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% BESTSPH
%
% DESCRIPTION
% This is equivalent to calling the BESTSPH macro function in CODE V.
%
% Surface	    Input: The desired surface (Integer)
% ZoomPos	    Input: The desired zoom position (Integer)
% MinHeight	    Input: Minimum y coordinate (omit center values of this limit) (Double)
% MaxHeight	    Input: Maximum y coordinate (Double)
% Curvature	    Output: The return value of the BESTSPH macro function.  It is the  
%               curvature of the best fitting sphere. (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    invoke(cv,'command',['res ' lens_to_load]) ;
    Surface = 3 ;
    ZoomPos = 1 ;
    MinHeight = 0 ;
    MaxHeight = 14.02598 ;

    clear Curvature ;
    [Curvature] = invoke(cv,'BESTSPH',Surface,ZoomPos,MinHeight,MaxHeight) ;

    disp(' ')
    disp(['Curvature is = ' num2str(Curvature)])
    disp(['Radius of Curvature = ' num2str((1/Curvature))])
    
    clear Surface ZoomPos MinHeight MaxHeight Curvature Radius ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ZERNIKE
%
% DESCRIPTION
% This is equivalent to calling the ZERNIKE macro function in CODE V.
%
% WaveNum	    Input: Number of the wavelength to use (not the value of the
%               wavelength) (Integer)
% FieldNum	    Input: Number of the field to use (not the value of the field) (Integer)
% ZoomPos	    Input: Zoom position to use. (Integer)
% CoefNum	    Input: Zernike coefficient to compute (1 to 37) (Integer)
% NumRays	    Input: Number of rays across the pupil diameter 
%               (0 uses a default of 21) (Integer)
% NumTerms	    Input: Number of Fringe Zernike terms to use in computing the   
%               coefficient. Zero uses a default of 37, the full Fringe 
%               Zernike set. (Integer)
% PupilType	    Input: An enum containing the pupil type. (Integer)
%                   0 ENP
%                   1 EXP
%                   2 EXS
% PolType	    Input: A number specifying if the polarization ray tracing is enabled
%               for this computation. See CODE V manual for explanation. (Integer)
% OutputType	Input: An enum containing the output type (Integer)
%                   0 Intensity
%                   1 Phase
% ZernType	    Input: An enum containing the Zernike type (Integer)
%                   0 ZFR
%                   1 ZRN
% Result	    Output: The return value is the Zernike coefficient; it is zero 
%               if a failure occurs. (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    invoke(cv,'command',['res ' lens_to_load]) ;
    invoke(cv,'command','thi si 0.01') ;
    
    WaveNum = 1 ;
    FieldNum = 1 ;
    ZoomPos = 1 ;
    CoefNum = 5 ;
    NumRays = 21 ;      % 0 is default of 21
    NumTerms = 28 ;     % Number of zernike terms to fit with
    PupilType = 0 ;      % (0) ENP, (1) EXP, (2) EXS
    PolType = 0 ;       % (0) polarization ray trace disabled, 
                        % see CODE V manual for other choices
    OutputType = 1 ;    % (0) Intensity (0),  Phase (1) zernike terms
    ZernType = 1 ;      % (0) ZFR or (1) ZRN

    clear result ;
    [result] = invoke(cv,'ZERNIKE',WaveNum,FieldNum,ZoomPos,CoefNum,...
                NumRays,NumTerms,PupilType,PolType,OutputType,ZernType) ;
    
    switch ZernType
        case 0
            ZernType = 'ZFR' ;
        otherwise    
            ZernType = 'ZRN' ;
    end
    
    disp(' ') 
    disp(['Zernike Coefficient Z' int2str(CoefNum) ' = ' num2str(result)...
                            ' For ' ZernType ' Type'])
            
    clear WaveNum FieldNum ZoomPos CoefNum NumRays NumTerms PupilType ...
                                                PolType OutputType ZernType result 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ZFRCOEF
%
% DESCRIPTION
% This is equivalent to calling the ZFRCOEF macro function in CODE V.
%
% WaveNum	    Input: Number of the wavelength to use (not the value of 
%               the wavelength) (Integer)
% FieldNum	    Input: Number of the field to use (not the value of the 
%               field) (Integer)
% ZoomPos	    Input: Zoom position to use. (Integer)
% CoefNum	    Input: Zernike coefficient to compute (1 to 37) (Integer)
% NumRays	    Input: Number of rays across the pupil diameter (0 uses a 
%               default of 21) (Integer)
% NumTerms	    Input: Number of Fringe Zernike terms to use in computing 
%               the coefficient.  Zero uses a default of 37, the 
%               full Fringe Zernike set.  (Integer)
% PupilType	    Input: An enum containing the pupil type.  (Integer)
%                   0 ENP
%                   1 EXP
%                   2 EXS
% Result	    Output: The return value is the Zernike coefficient; it is 
%               zero if a failure occurs. (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
    invoke(cv,'command',['res ' lens_to_load]) ;
    invoke(cv,'command','thi si 0.01') ;
    
    WaveNum = 1 ;
    FieldNum = 1 ;
    ZoomPos = 1 ;
    CoefNum = 4 ;
    NumRays = 21 ;      % 0 is default of 21
    NumTerms = 37 ;     % Fit uses default of 37 fringe terms
    PupilType = 0 ;      % (0) ENP, (1) EXP, (2) EXS

    clear result ;
    [result] = invoke(cv,'ZFRCOEF',WaveNum,FieldNum,ZoomPos,CoefNum,...
                                                        NumRays,NumTerms,PupilType) ;
    
    disp(' ')                                                
    disp(['Zernike Coefficient Z' int2str(CoefNum) ' = ' num2str(result)...
                                                ' For ZFR Type']) 
    
        clear WaveNum FieldNum ZoomPos CoefNum NumRays NumTerms PupilType result
                                              
                                                    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EVALZERN
%
% DESCRIPTION
% This is equivalent to calling the EVALZERN macro function in CODE V.
%
% WaveNum    	    Input: Number of the wavelength defined in ZERNIKE or 
%                   ZFRCOEF. (Integer)
% FieldNum	        Input: Number of the field point defined in ZERNIKE or 
%                   ZFRCOEF (Integer)
% ZoomPos	        Input: Zoom position defined in ZERNIKE or ZFRCOEF (Integer)
% X	                Input: X coordinate to be evaluated (Double)
% Y  	            Input: Y coordinate to be evaluated (Double)
% PolType	        Input: Number specifying whether polarization ray tracing  
%                   is enabled for this computation; matches the number defined 
%                   in ZERNIKE.  If you used ZFRCOEF, this must be 0. See CODE V 
%                   manual.(Integer)          
% OutputType	    Input: Enumeration of the output type.  Matches the type used in   
%                   ZERNIKE.  If you used ZFRCOEF, it must be phase. (Integer)
%                       0 Intensity
%                       1 Phase
% ZernType	        Input: Type of the Zernike polynomial.  Matches the expression 
%                   defined in ZERNIKE.  If you used ZFRCOEF, this must be 
%                   ZFR. (Integer)         
%                       0 ZFR
%                       1 ZRN
% Result	        Output: The value of the Zernike polynomial at the specified  
%                   coordinate. If the polynomial has not been defined with the  
%                   ZERNIKE or ZFRCOEFfunction, EVALZERN returns a value 
%                   of -1e10. (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
   invoke(cv,'command',['res ' lens_to_load]) ;
   invoke(cv,'command','thi si 0.01') ;
    
    ZernikeComputeChoice = 'ZFRCOEF' ;  % ZFRCOEF or ZERNIKE
    
    switch ZernikeComputeChoice

        case 'ZFRCOEF'
                WaveNum = 1 ;
                FieldNum = 1 ;
                ZoomPos = 1 ;
                CoefNum = 4 ;
                NumRays = 21 ;      % 0 is default of 21
                NumTerms = 37 ;     % Fit uses default of 37 fringe terms
                PupilType = 0 ;     % (0) ENP, (1) EXP, (2) EXS

                clear result ;
                [result] = invoke(cv,'ZFRCOEF',WaveNum,FieldNum,ZoomPos,CoefNum,...
                                                        NumRays,NumTerms,PupilType) ;
        
                X = 0.5 ;
                Y = 0.3 ;
                PolType = 0 ;       % Must be 0 polarization ray trace disabled
                OutputType = 1 ;    % Must be Phase (1) 
                ZernType = 0 ;      % Must be (0) ZFR

                clear result ;
                [result] = invoke(cv,'EVALZERN',WaveNum,FieldNum,ZoomPos,X,Y,...
                                                        PolType,OutputType,ZernType) ;
        case 'ZERNIKE'
            
                WaveNum = 1 ;
                FieldNum = 1 ;
                ZoomPos = 1 ;
                CoefNum = 5 ;
                NumRays = 21 ;      % 0 is default of 21
                NumTerms = 28 ;     % Number of zernike terms to fit with
                PupilType = 1 ;      % (0) ENP, (1) EXP, (2) EXS
                PolType = 0 ;       % (0) polarization ray trace disabled, 
                                    % see CODE V manual for other choices
                OutputType = 1 ;    % (0) Intensity (0),  Phase (1) zernike terms
                ZernType = 1 ;      % (0) ZFR or (1) ZRN

                clear result ;
                [result] = invoke(cv,'ZERNIKE',WaveNum,FieldNum,ZoomPos,CoefNum,...
                            NumRays,NumTerms,PupilType,PolType,OutputType,ZernType) ;
            
                X = 0.5 ;
                Y = 0.3 ;
                
                clear result ;
                [result] = invoke(cv,'EVALZERN',WaveNum,FieldNum,ZoomPos,X,Y,...
                                                        PolType,OutputType,ZernType) ;
                                                    
    end
    
    switch PupilType
        case 0
            Pupil = 'ENP' ;
        case 1
            Pupil = 'EXP' ;
        case 2
            Pupil = 'EXS' ;
        otherwise
            error('Incorrect Pupil Choice EVALZERN') 
    end
    
    
    disp(' ') 
    disp(['Zernike Z' int2str(CoefNum) ' = ' num2str(result)...
         ' for ' ZernikeComputeChoice '. Evaluation at (X,Y) = '...
             '(' num2str(X) ',' num2str(Y) ')' ' in the ' Pupil ])      
                                        
    clear WaveNum FieldNum ZoomPos CoefNum NumRays NumTerms PupilType ...
            X Y PolType OutputType ZernType result Pupil ZernikeComputeChoice 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FITERROR
%
% DESCRIPTION
% This is equivalent to calling the FITERROR macro function in CODE V.
%
% WaveNum    	Input: Number of the wavelength defined in ZERNIKE or ZFRCOEF (Integer)
% FieldNum	    Input: Number of the field point defined in ZERNIKE or ZFRCOEF (Integer)
% ZoomPos	    Input: Zoom position defined in ZERNIKE or ZFRCOEF (Integer)
% PolType	    Input: Number specifying whether polarization ray tracing is enabled  
%               for this computation; matches the number defined in ZERNIKE. If you  
%               used ZFRCOEF, this must be 0. See CODE V manual.(Integer)
% OutputType	Input: Enumeration of the output type.  Matches the type used   
%               in ZERNIKE.  If you used ZFRCOEF, it must be phase.(Integer)
%                       0 Intensity
%                       1 Phase
% ZernType	    Input: Type of the Zernike polynomial.  Matches the expression defined   
%               in ZERNIKE. If you used ZFRCOEF, this must be ZFR. (Integer)
%                       0 ZFR
%                       1 ZRN
% Result	    Output: The RMS fit error of the Zernike polynomial.  If the  
%               polynomial has not been defined with the ZERNIKE or ZFRCOEF function, 
%               EVALZERN returns a value of -1. (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
   invoke(cv,'command',['res ' lens_to_load]) ;
   invoke(cv,'command','thi si 0.01') ;
    
    ZernikeComputeChoice = 'ZERNIKE' ;  % ZFRCOEF or ZERNIKE
    
    switch ZernikeComputeChoice

        case 'ZFRCOEF'
                WaveNum = 1 ;
                FieldNum = 1 ;
                ZoomPos = 1 ;
                CoefNum = 4 ;
                NumRays = 21 ;      % 0 is default of 21
                NumTerms = 28 ;     % Fit uses default of 37 fringe terms
                PupilType = 1 ;     % (0) ENP, (1) EXP, (2) EXS

                clear result ;
                [result] = invoke(cv,'ZFRCOEF',WaveNum,FieldNum,ZoomPos,CoefNum,...
                                                        NumRays,NumTerms,PupilType) ;
        
                PolType = 0 ;       % Must be 0 polarization ray trace disabled
                OutputType = 1 ;    % Must be Phase (1) 
                ZernType = 0 ;      % Must be (0) ZFR

                clear result ;
                [result] = invoke(cv,'FITERROR',WaveNum,FieldNum,ZoomPos,...
                                                        PolType,OutputType,ZernType) ;
        case 'ZERNIKE'
            
                WaveNum = 1 ;
                FieldNum = 1 ;
                ZoomPos = 1 ;
                CoefNum = 5 ;
                NumRays = 21 ;      % 0 is default of 21
                NumTerms = 28 ;     % Number of zernike terms to fit with
                PupilType = 1 ;      % (0) ENP, (1) EXP, (2) EXS
                PolType = 0 ;       % (0) polarization ray trace disabled, 
                                    % see CODE V manual for other choices
                OutputType = 1 ;    % (0) Intensity (0),  Phase (1) zernike terms
                ZernType = 0 ;      % (0) ZFR or (1) ZRN

                clear result ;
                [result] = invoke(cv,'ZERNIKE',WaveNum,FieldNum,ZoomPos,CoefNum,...
                            NumRays,NumTerms,PupilType,PolType,OutputType,ZernType) ;
            
                clear result ;        
                [result] = invoke(cv,'FITERROR',WaveNum,FieldNum,ZoomPos,...
                                                        PolType,OutputType,ZernType) ;
        otherwise
                error('Invalide FITERROR CHOICE')
                                                  
    end
    
    disp(' ') 
    disp(['Zernike Fit Error = ' num2str(result)]) 
     
    clear WaveNum FieldNum ZoomPos CoefNum NumRays NumTerms PupilType ...
                                           X Y PolType OutputType ZernType result 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% NORMRADIUS
%
% DESCRIPTION
% This is equivalent to calling the NORMRADIUS macro function in CODE V.
%
% WaveNum	    Input: Number of the wavelength defined in ZERNIKE or 
%               ZFRCOEF (Integer)
% ZoomPos	    Input: Zoom position defined in ZERNIKE or ZFRCOEF (Integer)
% PolType	    Input: Number specifying whether polarization ray tracing is enabled 
%               for this computation; matches the number defined in ZERNIKE.  If 
%               you used ZFRCOEF, this must be 0. See CODE V manual. (Integer)
% OutputType	Input: Enumeration of the output type.  Matches the type used 
%               in ZERNIKE. If you used ZFRCOEF, it must be phase.
%                       0 Intensity
%                       1 Phase
% ZernType	    Input: Type of the Zernike polynomial.  Matches the expression    
%               defined in ZERNIKE. If you used ZFRCOEF, this must be 
%               ZFR. (Integer)
%                       0 ZFR
%                       1 ZRN
% Result	    Output: The normalization radius of the Zernike polynomial.  If the  
%               polynomial has not been defined with the ZERNIKE or ZFRCOEF 
%               function, EVALZERN returns a value of -1.  (Double)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
   invoke(cv,'command',['res ' lens_to_load]) ;
   invoke(cv,'command','thi si 0.01') ;
    
    ZernikeComputeChoice = 'ZFRCOEF' ;  % ZFRCOEF or ZERNIKE
    
    switch ZernikeComputeChoice

        case 'ZFRCOEF'
                WaveNum = 1 ;
                FieldNum = 1 ;
                ZoomPos = 1 ;
                CoefNum = 4 ;
                NumRays = 21 ;      % 0 is default of 21
                NumTerms = 28 ;     % Fit uses default of 37 fringe terms
                PupilType = 0 ;     % (0) ENP, (1) EXP, (2) EXS

                clear result ;
                [result] = invoke(cv,'ZFRCOEF',WaveNum,FieldNum,ZoomPos,CoefNum,...
                                                        NumRays,NumTerms,PupilType) ;
        
                PolType = 0 ;       % Must be 0 polarization ray trace disabled
                OutputType = 1 ;    % Must be Phase (1) 
                ZernType = 0 ;      % Must be (0) ZFR

                clear result ;
                [result] = invoke(cv,'NORMRADIUS',WaveNum,FieldNum,ZoomPos,...
                                                        PolType,OutputType,ZernType);
                                                    
        case 'ZERNIKE'
            
                WaveNum = 1 ;
                FieldNum = 1 ;
                ZoomPos = 1 ;
                CoefNum = 5 ;
                NumRays = 21 ;      % 0 is default of 21
                NumTerms = 28 ;     % Number of zernike terms to fit with
                PupilType = 0 ;      % (0) ENP, (1) EXP, (2) EXS
                PolType = 0 ;       % (0) polarization ray trace disabled, 
                                    % see CODE V manual for other choices
                OutputType = 1 ;    % (0) Intensity (0),  Phase (1) zernike terms
                ZernType = 0 ;      % (0) ZFR or (1) ZRN

                clear result ;
                [result] = invoke(cv,'ZERNIKE',WaveNum,FieldNum,ZoomPos,CoefNum,...
                            NumRays,NumTerms,PupilType,PolType,OutputType,ZernType) ;
            
                clear result ;
                [result] = invoke(cv,'NORMRADIUS',WaveNum,FieldNum,ZoomPos,...
                                                        PolType,OutputType,ZernType);
       otherwise
                error('Invalide NORMRADIUS CHOICE')
                                                  
    end
    
    switch ZernikeComputeChoice
        case 'ZFRCOEF'
            z = 'ZFR' ;
        case 'ZERNIKE'
            z = 'ZRN' ;
    end
    
    disp(' ')       
    disp(['Normalization Radius = ' num2str(result) ' for ' z ...
                                       ' Z' int2str(CoefNum)])
    
    clear WaveNum FieldNum ZoomPos CoefNum NumRays NumTerms PupilType ...
                            z ZernikeComputeChoice X Y PolType OutputType ZernType result 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% RAYRSI
%
% DESCRIPTION
% This is equivalent to calling the RAYRSI macro function in CODE V.
%
% ZoomPos	    Input: Zoom position to ray trace.  If the value is < 1 or >( num z),  
%               then zoom position 1 is used. (Integer)
% WaveNum	    Input: The wavelength number to use.  If the value is < 1 or > (num w), 
%               then the reference wavelength is used. (Integer)
% FieldNum	    Input: The field number to use.  If the value > 0 then the last two  
%               elements of input are ignored and the field number is used.  If  
%               the value is < 1 or > (num f), then it is ignored and input(3) and 
%               input(4) are used. (Integer)
% RefSurface	Input: The reference surface. (Integer)  
% Input	        Input: A 4 element input vector containing the relative pupil and 
%               relative field coordinates. (Double)
% Result	    Output: The return value is the surface number at which the ray trace   
%               failed.  A return of 0 indicates success.  RER database access -- so a
%               catastrophic failure is needed, i.e. missed surface) (Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
   invoke(cv,'command',['res ' lens_to_load]) ;

   ZoomPos = 1 ;
   WaveNum = 2 ;
   FieldNum = 3 ;
   RefSurface = 1 ;
   Input = zeros(4,1) ;
   Input(1) = 1 ;   % x either coordinate on surface or relative pupil (1=full pupil)
   Input(2) = 0 ;   % y either coordinate on surface or relative pupil
   Input(3) = 1 ;   % x relative field (1=the full of that field)
   Input(4) = 1 ;   % y relative field 
  
   clear result ;
   [result] = invoke(cv,'RAYRSI',ZoomPos,WaveNum,FieldNum,RefSurface,Input) ;

   clear C
   C(1,1) = {'x-coordinate'} ;
   C(2,1) = {'y-coordinate'} ;
   C(3,1) = {'x relative field'} ;
   C(4,1) = {'y relative field'} ;
   
   for ii = 1:4
       C(ii,2) = {Input(ii)} ;
   end
   
   disp(' ')
   disp('RAYRSI Trace')
   disp(C)
   disp(['RAYRSI result = ' int2str(result)])
   
   clear result ZoomPos WaveNum FieldNum RefSurface Input ii C;
end   
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% RAYSIN
%
% DESCRIPTION
% This is equivalent to calling the RAYSIN macro function in CODE V.
%
% ZoomPos	     Input: Zoom position to ray trace.  If the value is < 1 or >( num z), 
%                then zoom position 1 is used. (Integer)
% WaveNum	     Input: The wavelength number to use.  If the value is < 1 or >  
%                (num w), then the reference wavelength is used. (Integer)
% X              Input: The X ray coordinates on the tangent plane for 
%                surface 1. (Double)
% Y  	         Input: The Y ray coordinates on the tangent plane for 
%                surface 1. (Double)
% XDirTan        Input: Input: The X direction tangents in object space. (Double)
% YDirTan	     Input: The Y direction tangents in object space. (Double)
% Result	     Output: The return value is the surface number at which the ray trace 
%                failed.  A return of 0 indicates success. (RER database check
%                -- so a catastrophic failure is needed, i.e. missed
%                surface) (Integer)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
   invoke(cv,'command',['res ' lens_to_load]) ;
   invoke(cv,'command','ins s1') ;

   ZoomPos = 1 ;   
   WaveNum = 2 ;   
   X = 16.6667 ;
   Y = 0.0 ;
   XDirTan = 0 ;
   YDirTan = 0 ;
   clear result ;
   [result] = invoke(cv,'RAYSIN',ZoomPos,WaveNum,X,Y,XDirTan,YDirTan) ;

   clear C
   C(1,1) = {'x-coordinate'} ;
   C(2,1) = {'y-coordinate'} ;
   C(3,1) = {'x direction tan'} ;
   C(4,1) = {'y direction tan'} ;
   
   C(1,2) = {X} ;
   C(2,2) = {Y} ;
   C(3,2) = {XDirTan} ;
   C(4,2) = {YDirTan} ;
   
   disp(' ')
   disp('RAYSIN Trace')
   disp('Input Vector At Tangent Plane To Surface 1') 
   disp(C)
   disp(['RAYSIN result = ' int2str(result)])
   
   clear result ZoomPos WaveNum X Y XDirTan YDirTan C ;
   
end   
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAYTRA
% 
% DESCRIPTION
% This is equivalent to calling the RAYTRA macro function in CODE V.
%
% ZoomPos	Input: Zoom position to ray trace.  If the value is < 1 or >( num z), then 
%           zoom position 1 is used. (Integer)
% WaveNum	Input: The wavelength number to use.  If the value is < 1 or > (num w), then 
%           the reference wavelength is used. (Integer)
% AperCheck	Input: Flag to check apertures.  0 means do not check apertures, non-zero 
%           means to check apertures. (Integer)
% Input	    Input: A 4 element vector containing the coordinates on the first tangent  
%           plane and the direction tangents in object space. (Double)
% Output	Output: An 8 element array containing the image surface ray coordinates  
%           (X,Y,Z),the direction cosines (L,M,N), optical path from surface 1 tangent 
%           plane to the image surface, and transmission of the exiting array.
%           (Double)
% Result	Output: The return value is as follows: (Integer)
%               0 Succesful trace
%               -n TIR or diffraction into evanescent order at surface n
%                n Ray blocked at surface n
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (1)
    
   invoke(cv,'command',['res ' lens_to_load]) ;
   invoke(cv,'command','ins s1') ;
   
   ZoomPos = 1 ;
   WaveNum = 2 ;
   AperCheck = 1 ;
   Input = zeros(4,1) ;
   Input(1) = 16.5 ;    % x position on surface 1 tangent plane
   Input(2) = 0 ;       % y position on surface 1 tangent plane
   Input(3) = 0.0 ;     % x direction tangent
   Input(4) = 0.1 ;     % y direction tangent 
   Output = zeros(8,1) ;
  
   clear result ;
   [result,Output] = invoke(cv,'RAYTRA',ZoomPos,WaveNum,AperCheck,Input,Output) ;

   Inputcell(1,1) = {'Zoom Position'} ;
   Inputcell(2,1) = {'Wavelength'} ;
   Inputcell(3,1) = {'Apercheck'} ;
   
   Inputcell(1,2) = {ZoomPos} ;
   Inputcell(2,2) = {WaveNum} ;
   Inputcell(3,2) = {AperCheck} ;
   
   Inputcell(4,1) = {'x'} ;
   Inputcell(5,1) = {'y'} ;
   Inputcell(6,1) = {'x tangent'} ;
   Inputcell(7,1) = {'y tangent'} ;
   
   for ii = 1:4
       Inputcell(ii+3,2) = {Input(ii)} ;
   end
   
   Outputcell(1,1) = {'x'} ;
   Outputcell(2,1) = {'y'} ;
   Outputcell(3,1) = {'z'} ;
   Outputcell(4,1) = {'l'} ;
   Outputcell(5,1) = {'m'} ;
   Outputcell(6,1) = {'n'} ;
   Outputcell(7,1) = {'Optical Path'} ;
   Outputcell(8,1) = {'Transmission'} ; 
   
   for ii = 1:8
       Outputcell(ii,2) = {Output(ii)} ;
   end
   
   disp(' ')
   disp('Input Vector At Tangent Plane To Surface 1') 
   disp(Inputcell)
   disp('Image Plane Coordinates')
   disp(Outputcell)
   disp(['RAYSIN result = ' int2str(result)])
   
   clear result Output ZoomPos WaveNum AperCheck Input Inputcell Outputcell ii;
    
end   
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% StopCodeV
%
% DESCRIPTION
% StopCodeV must be called when you are done running the session of CODE V.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

invoke(cv,'StopCodeV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cv paths lens_to_load ;;
