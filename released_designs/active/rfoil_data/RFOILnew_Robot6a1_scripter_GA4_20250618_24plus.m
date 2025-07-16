%RFoil Airfoil Script Robot
%This script sends key-press commands to RFOIL using JAVA Robot
%If you click on a word processor, you can see the commands through the
%first alfa case

%v4b (10/23/2024 dcm): Adjusted code to accept a table input to run
%multiple airfoils for the fitted FFA-W3's.  
%4c: (10/23/24): adding ability to handle rough condition in vpar menu.
%5a (12-3-2024) %Added a temp output file switch.  Have RFoil output as a 'temp_out.dat' file
%               that you know will work, and then change the file name to what you want in matlab.
%TODO: 
% -add flag to make the initialization AOA cases skipped or shorter.
% -make sure it works without the user click over...


%steps to use:
%Pre-run:  set AOA Low, High, delta; set Re# list; set ...
%AirfoilPolarOutFile; set airfoil name

% 1.)  Start RFOIL
% 2.) Run this script
% 3.) Click on RFOIL window within 'UserClickOverDelay' seconds.  
% DON'T TOUCH ANYTHING UNTIL THE SCRIPT IS FINISHED
%4.)  Learn from all of the data!



clear

TimeStartDays = now; %'now' gives the current time as the serial number of days since the beginning of this century

USER_CLICKOVER_INPUT = 0; %This determines if the script waits for the user to click on the matlab script, press ENTER, then click on the RFoil window
%or if it just waits for the AOASEQ_DELAY
AOASEQ_DELAY = 40; %amount of delay to allow AOA loop to run, most runs are about 20 seconds, this should be very conservative

RFoil_Version = 'New';%3/31/2025 RFoil version, New = 3.04, Old = 1.0


%% User Inputs

% #     # tau,  CL,  spn,     Re
% #     [0.15, 1.5, 1.00, 10.0e6, ],
% #     [0.18, 1.5, 1.00, 10.0e6, ],
% #     [0.21, 1.5, 1.00, 12.0e6, ],
% #     [0.24, 1.4, 0.85, 13.0e6, ],
% #     [0.27, 1.3, 0.55, 16.0e6, ],
% #     [0.30, 1.2, 0.50, 18.0e6, ],
% #     [0.33, 1.2, 0.35, 16.0e6, ],
% #     [0.36, 1.2, 0.20, 13.0e6, ],


%Table(Names) = Case#, AirfoilPolarOutFile, AirfoilCoordFile, AirfoilName, Clean or Rough (C or R), AOA_Low, AOA_High, AOA_delta, Mach, DRAG_NewModel, [ReNum1, ReNum2, etc]
AirfoilCell =  ...
{ 1, 'final_airfoil_24_r13_cleanPol_RfoilNew6a_DragOff.dat', 'final_airfoil_24.dat', 'final_airfoil_24', 'C', -24, 20, 1, 0.0, 'Off', [13e6];
  2, 'final_airfoil_24_r13_roughPol_RfoilNew6a_DragOff.dat', 'final_airfoil_24.dat', 'final_airfoil_24', 'R', -24, 20, 1, 0.0, 'Off', [13e6];
  3, 'final_airfoil_24_r13_cleanPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_24.dat', 'final_airfoil_24', 'C', -24, 20, 1, 0.0, 'On' , [13e6];
  4, 'final_airfoil_24_r13_roughPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_24.dat', 'final_airfoil_24', 'R', -24, 20, 1, 0.0, 'On' , [13e6];
     5, 'final_airfoil_27_r16_cleanPol_RfoilNew6a_DragOff.dat', 'final_airfoil_27.dat', 'final_airfoil_27', 'C', -24, 20, 1, 0.0, 'Off', [16e6];
  6, 'final_airfoil_27_r16_roughPol_RfoilNew6a_DragOff.dat', 'final_airfoil_27.dat', 'final_airfoil_27', 'R', -24, 20, 1, 0.0, 'Off', [16e6];
  7, 'final_airfoil_27_r16_cleanPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_27.dat', 'final_airfoil_27', 'C', -24, 20, 1, 0.0, 'On' , [16e6];
  8, 'final_airfoil_27_r16_roughPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_27.dat', 'final_airfoil_27', 'R', -24, 20, 1, 0.0, 'On' , [16e6];
     9, 'final_airfoil_30_r18_cleanPol_RfoilNew6a_DragOff.dat', 'final_airfoil_30.dat', 'final_airfoil_30', 'C', -24, 20, 1, 0.0, 'Off', [18e6];
  10, 'final_airfoil_30_r18_roughPol_RfoilNew6a_DragOff.dat', 'final_airfoil_30.dat', 'final_airfoil_30', 'R', -24, 20, 1, 0.0, 'Off', [18e6];
  11, 'final_airfoil_30_r18_cleanPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_30.dat', 'final_airfoil_30', 'C', -24, 20, 1, 0.0, 'On' , [18e6];
  12, 'final_airfoil_30_r18_roughPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_30.dat', 'final_airfoil_30', 'R', -24, 20, 1, 0.0, 'On' , [18e6];
     13, 'final_airfoil_33_r16_cleanPol_RfoilNew6a_DragOff.dat', 'final_airfoil_33.dat', 'final_airfoil_33', 'C', -24, 20, 1, 0.0, 'Off', [16e6];
  14, 'final_airfoil_33_r16_roughPol_RfoilNew6a_DragOff.dat', 'final_airfoil_33.dat', 'final_airfoil_33', 'R', -24, 20, 1, 0.0, 'Off', [16e6];
  15, 'final_airfoil_33_r16_cleanPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_33.dat', 'final_airfoil_33', 'C', -24, 20, 1, 0.0, 'On' , [16e6];
  16, 'final_airfoil_33_r16_roughPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_33.dat', 'final_airfoil_33', 'R', -24, 20, 1, 0.0, 'On' , [16e6];
     17, 'final_airfoil_36_r13_cleanPol_RfoilNew6a_DragOff.dat', 'final_airfoil_36.dat', 'final_airfoil_36', 'C', -24, 20, 1, 0.0, 'Off', [13e6];
  18, 'final_airfoil_36_r13_roughPol_RfoilNew6a_DragOff.dat', 'final_airfoil_36.dat', 'final_airfoil_36', 'R', -24, 20, 1, 0.0, 'Off', [13e6];
  19, 'final_airfoil_36_r13_cleanPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_36.dat', 'final_airfoil_36', 'C', -24, 20, 1, 0.0, 'On' , [13e6];
  20, 'final_airfoil_36_r13_roughPol_RfoilNew6a_DragOn.dat' , 'final_airfoil_36.dat', 'final_airfoil_36', 'R', -24, 20, 1, 0.0, 'On' , [13e6];}
   % 3, 'doe30_r6_cleanPol.dat', 'jan-9-doe_30.dat', 'jan-9-doe_30', 'C', -5, 0, 1, 0.0, [6e6];
   % 4, 'doe30_r6_roughPol.dat', 'jan-9-doe_30.dat', 'jan-9-doe_30', 'R', -5, 0, 1, 0.0, [6e6];
   % 5, 'doe30_r10_cleanPol.dat', 'jan-9-doe_30.dat', 'jan-9-doe_30', 'C', -5, 0, 1, 0.0, [10e6];
   % 6, 'doe30_r10_roughPol.dat', 'jan-9-doe_30.dat', 'jan-9-doe_30', 'R', -5, 0, 1, 0.0, [10e6];}
%Seems like output file is limited to 12 characters in RFoil2007??

RFoil_Robot(AirfoilCell, USER_CLICKOVER_INPUT, AOASEQ_DELAY, RFoil_Version);

function RFoil_Robot(AFCell, USER_CLICKOVER_INPUT, AOASEQ_DELAY, RFoil_Version)

    FileNameOutTemp = 'TempOut.dat';%This is always the file RFoil will output, and then file name gets changed to the desired output file.

    for i = 1:1:size(AFCell,1)
        airfoilNum =    int16(AFCell{i,1});%must be integer
        File_AirfoilOut = AFCell{i,2};
        File_AirfoilIn = AFCell{i,3};
        Airfoil_Name = AFCell{i,4};
        RoughCleanFlag = AFCell{i,5};
        AOA_Low = AFCell{i,6};
        AOA_High = AFCell{i,7};
        AOA_Delta = AFCell{i,8};
        Mach = AFCell{i,9};
        DRAG_NewModel = AFCell{i,10};
        ReNumList = AFCell{i,11};
        
    
        if(RFoil_Version == 'Old')
            RFoil_Robot_innerloop(airfoilNum, FileNameOutTemp, File_AirfoilIn, Airfoil_Name,RoughCleanFlag,AOA_Low, AOA_High, AOA_Delta, Mach, ReNumList, USER_CLICKOVER_INPUT, AOASEQ_DELAY )
        elseif(RFoil_Version == 'New')
            RFoilnew_Robot_innerloop(airfoilNum, FileNameOutTemp, File_AirfoilIn, Airfoil_Name,RoughCleanFlag,AOA_Low, AOA_High, AOA_Delta, Mach, ReNumList, DRAG_NewModel)
        end
        copyfile(FileNameOutTemp, File_AirfoilOut);
        delete(FileNameOutTemp);

    end
end

%%RFoil function for new RFoil v3+
function RFoilnew_Robot_innerloop(airfoilNum, File_AirfoilOut, File_AirfoilIn, Airfoil_Name,RoughCleanFlag,AOA_Low, AOA_High, AOA_Delta, Mach, ReNumList, DRAG_NewModel)
  disp(strcat('Case airfoilNum: ',num2str(airfoilNum)))
%Run Flags:
DEBUG = 0;%RFoilnew v6a, debug not used yet, used to add pause after each major step
AOA_INIT_EACH_RENUM = 0;%Sets if extra AoA initialization steps are needed for each ReNum.  Not needed if AOA_Low is > -5 deg. or if not doign multiple ReNum
                    %RFoilnew v6a, not used yet


%DRAG_NewModel = 'OFF'; %New Rfoil Drag model 'OFF' or 'ON'

%File_AirfoilIn = '6_to_1_test_airfoil_24.dat'; %coordinates file
%Airfoil_Name = 'doe_24_6to1';
%AOA_Delta = 1.0;
%Mach = 0.0;
%ReNumList(1) = 13e6;
%USER_CLICKOVER_INPUT, AOASEQ_DELAY; not used with new RFoil version!

 %AOA_Low = -10.0;%Low AoA for aseq
 %AOA_High = 22.0;%High AoA for aseq

%filename_in = 'rfoil_clean_DragOff_TestIn.inp'

%filename_out_RFoilInp = 'RFoilNew_TestIn.inp'

outputFilePath_StdOut = 'std_out.txt';

N_factor_Clean = 9.0;%n=9 is typically 'clean, 0.070 Turbulence)
    XTR_S_Clean = 1.0;%XTR/c, S = suction surface (upper) transition location, 1.0 for clean. 0.0 for Leading Edge
    XTR_P_Clean = 1.0;%XTR/c, P = pressure Surface (lower)
    N_factor_Rough = 3.0;%n=9 is typically 'clean, 0.070 Turbulence)
    XTR_S_Rough = 0.05;%S = suction surface (upper)
    XTR_P_Rough = 0.05;%P = pressure Surface (lower)

   ChangePanel = 1;%if 0, then paneling will be default (160), if 1 then paneling set by PanelNum
   PanelNum = 240;%Note, if you increase this above 160, should increase AOASEQ_DELAY > 30s

   % AirfoilPolarOutFile = [File_AirfoilOut];
   % AirfoilCoordFile = File_AirfoilIn;
   % AirfoilName = Airfoil_Name;

   %AirfoilPolarOutFile = [File_AirfoilOut];
   AirfoilCoordFile = File_AirfoilIn;
   AirfoilName = Airfoil_Name;

   AOA_delta = AOA_Delta;

MachInit = Mach;
ReInit = ReNumList(1);

%dcm 10-23-2024 adding roughness/clean case
if RoughCleanFlag == 'C' %Clean case
    N_Factor_VPAR = N_factor_Clean;%n=9 is typically 'clean, 0.075 Turbulence)
    XTR_S = XTR_S_Clean;%S = suction surface (upper)
    XTR_P = XTR_P_Clean;%P = pressure Surface (lower)
elseif RoughCleanFlag == 'R' %Rough Case
    N_Factor_VPAR = N_factor_Rough;%n=9 is typically 'clean, 0.075 Turbulence)
    XTR_S = XTR_S_Rough;%S = suction surface (upper)
    XTR_P = XTR_P_Rough;%P = pressure Surface (lower)
else 
    disp('WARNING: RoughCleanFlag not equal to C or R in input table. Defaulting to Clean.')
    N_Factor_VPAR = N_factor_Clean;%n=9 is typically 'clean, 0.075 Turbulence)
    XTR_S = XTR_S_Clean;%XTR/c, S = suction surface (upper) transition location, 1.0 for clean. 0.0 for Leading Edge
    XTR_P = XTR_P_Clean;%XTR/c, P = pressure Surface (lower)
pause; 
end;

%Open File
FileOut_Rfoil = File_AirfoilOut;
%FileOut_Cln_Rfoil = 'TempOut_Cln3.dat';
%FileOut_Rgh_Rfoil = 'TempOut_Rgh3.dat';

if (exist(FileOut_Rfoil));
delete(FileOut_Rfoil);
end
% if (exist(FileOut_Rgh_Rfoil));
% delete(FileOut_Rgh_Rfoil);
% end

clear linesCellArray2;
linesCellArray2{1,1} =       'VERs'                           ;
linesCellArray2{end+1} =       strcat('Drag'," ",DRAG_NewModel) ;
linesCellArray2{end+1} =     'BACK'                           ;
linesCellArray2{end+1} =     'load'                           ;
linesCellArray2{end+1} =       AirfoilCoordFile   ;
linesCellArray2{end+1} =         Airfoil_Name                ;
linesCellArray2{end+1} =     'ppar'                           ;
linesCellArray2{end+1} =     'n'                              ;
linesCellArray2{end+1} =       num2str(PanelNum);
linesCellArray2{end+1} =     '';                         
linesCellArray2{end+1} =     'oper'                           ;
linesCellArray2{end+1} =     'visc'                           ;
linesCellArray2{end+1} =    num2str(ReInit)                        ;
linesCellArray2{end+1} =    'Mach'                           ;
linesCellArray2{end+1} =    num2str(MachInit)                              ;
linesCellArray2{end+1} =    'vpar'                           ;
linesCellArray2{end+1} =    'N'                              ;
linesCellArray2{end+1} =    num2str(N_Factor_VPAR);
  linesCellArray2{end+1} =  'xtr'                            ;
  linesCellArray2{end+1} =  num2str(XTR_S)                          ;
  linesCellArray2{end+1} =  num2str(XTR_P)                           ;
 linesCellArray2{end+1} =       '';                      
 linesCellArray2{end+1} =   'pacc';
 %linesCellArray2{end+1} =   'doe24_r12_cln_RfoilNew4TEST.dat';
  linesCellArray2{end+1} =   FileOut_Rfoil;%temp file out
 linesCellArray2{end+1} =     '';                        
 linesCellArray2{end+1} =   'aseq';
 linesCellArray2{end+1} =   num2str(AOA_Low);
 linesCellArray2{end+1} =   num2str(AOA_High);
 linesCellArray2{end+1} =   num2str(AOA_delta);
linesCellArray2{end+1} =    'pacc'                           ;%end pacc file store
 %Roug case
% linesCellArray2{end+1} =   'vpar'                           ;%prepare for rough case
%  linesCellArray2{end+1} =   'N'                              ;
%  linesCellArray2{end+1} =   num2str(N_factor_Rough);
%   linesCellArray2{end+1} =  'xtr'                            ;
%   linesCellArray2{end+1} =  num2str(XTR_S_Rough)                          ;
%   linesCellArray2{end+1} =  num2str(XTR_P_Rough)                           ;
%   linesCellArray2{end+1} =     '';                       
%   linesCellArray2{end+1} =  'pacc'                           ;
%   %linesCellArray2{end+1} =  'doe24_r12_rgh_RfoilNew4TEST.dat';
%     linesCellArray2{end+1} =  FileOut_Rgh_Rfoil;
%   linesCellArray2{end+1} =    '';                        
%   linesCellArray2{end+1} =  'aseq'                           ;
%  linesCellArray2{end+1} =   num2str(AOA_Low);
%  linesCellArray2{end+1} =   num2str(AOA_High);
%  linesCellArray2{end+1} =   num2str(AOA_delta);
%   linesCellArray2{end+1} =  'pacc'                           ;%end rough pacc file store
  linesCellArray2{end+1} =    '';                        %go up one level to quit
  linesCellArray2{end+1} =  'quit'                           ;

%MultipleTextLinesToFile(linesCellArray,filename_out_RFoilInp);
MultipleTextLinesToFile(linesCellArray2,'Temp_RFoilNew_TestIn3.inp');

%Run RFoil
%runExecutable('RFoil.exe < RFoilNew_TestIn.inp')
%runExecutableWithArgs('RFoil.exe', '"< RFoilNew_TestIn.inp"', '"argument with spaces"', 'arg3');

runExecutableWithArgs('RFoil.exe',outputFilePath_StdOut, '< Temp_RFoilNew_TestIn3.inp');

%End Main code
end

%% Functions

function MultipleTextLinesToFile(linesToWrite,filename)
% Define multiple lines of text in a cell array
% linesToWrite = {
%     'This is the first line of text.',
%     'This is the second line of text.',
%     'This is the third line of text.',
%     'This is the fourth line of text.'
% };

% Specify the filename
%filename = 'output_multiple_lines.txt';

% Open the file for writing (this will create the file if it does not exist)
fileID = fopen(filename, 'w');

% Check if the file opened successfully
if fileID == -1
    error('Failed to open the file for writing.');
end

% Write each line to the file
for i = 1:length(linesToWrite)
    fprintf(fileID, '%s\n', linesToWrite{i});
end

% Close the file
fclose(fileID);

disp('Multiple lines of text have been written to the file successfully.');
end

function runExecutableWithArgs(executablePath, outputFilePath, varargin)
    % runExecutableWithArgs Runs a specified Windows executable with additional arguments.
    %   runExecutableWithArgs(executablePath, arg1, arg2, ...) runs the executable
    %   located at executablePath with the specified additional command line arguments.
    %use: runExecutableWithArgs('C:\Path\To\YourExecutable.exe', 'arg1', 'arg2', 'arg3');
        %If any of the arguments contain spaces, you may need to wrap them in quotes when calling the function. For example:
        %runExecutableWithArgs('C:\Path\To\YourExecutable.exe', 'arg1', '"argument with spaces"', 'arg3');


    % Check if the executable exists
    if exist(executablePath, 'file') ~= 2
        error('The specified executable does not exist: %s', executablePath);
    end

    % Construct the command string
    %command = ['"' executablePath '"'];
command = [executablePath];

    % Append additional arguments
    for i = 1:length(varargin)
        command = [command ' ' varargin{i}];
    end

    % Run the executable with the constructed command
    disp('Starting RFoil in background')
    [status, cmdout] = system(command);

% Save the command line output to a file
    fileID = fopen(outputFilePath, 'w');
    if fileID == -1
        error('Failed to open output file: %s', outputFilePath);
    end
    fprintf(fileID, '%s', cmdout);
    fclose(fileID);

    % Check the status of the command
    if status == 0
        disp('Executable ran successfully.');
    else
        disp('Error running executable:');
        disp(['Status: ' num2str(status)]);
        disp(cmdout);
    end
end
    

%% RFoil old robot

function RFoil_Robot_innerloop(airfoilNum, File_AirfoilOut, File_AirfoilIn, Airfoil_Name,RoughCleanFlag,AOA_Low, AOA_High, AOA_Delta, Mach, ReNumList, USER_CLICKOVER_INPUT, AOASEQ_DELAY )
  disp(strcat('Case airfoilNum: ',num2str(airfoilNum)))
%Run Flags:
DEBUG = 0;
AOA_INIT_EACH_RENUM = 0;%Sets if extra AoA initialization steps are needed for each ReNum.  Not needed if AOA_Low is > -5 deg. or if not doign multiple ReNum.

    N_factor_Clean = 9.0;%n=9 is typically 'clean, 0.070 Turbulence)
    XTR_S_Clean = 1.0;%XTR/c, S = suction surface (upper) transition location, 1.0 for clean. 0.0 for Leading Edge
    XTR_P_Clean = 1.0;%XTR/c, P = pressure Surface (lower)
    N_factor_Rough = 3.0;%n=9 is typically 'clean, 0.070 Turbulence)
    XTR_S_Rough = 0.05;%S = suction surface (upper)
    XTR_P_Rough = 0.05;%P = pressure Surface (lower)

    ChangePanel = 0;%if 0, then paneling will be default (160), if 1 then paneling set by PanelNum
   PanelNum = 200;%Note, if you increase this above 160, should increase AOASEQ_DELAY > 30s

AirfoilPolarOutFile = [File_AirfoilOut]
%final output will be 'airfoilOut#_r010.dat'
%where 'r010' is the Reynolds number/10^5, so the above would be Re  1e6
%the '#' is an integer.  The number is cycled if an existing file is found.

%Airfoil file to load and name to store it as
    AirfoilCoordFile = File_AirfoilIn;
    AirfoilName = Airfoil_Name;

%The operating cl range (between negative and positive cl_max) tends
%to be a function of Reynolds number, with the range increasing as RENum
%increases.  Therefore, two ranges of angle of attack are allowed.
%More can easily be added as needed.  
%The low_RE cases are anything less than a Reynodls number of AOA_ReNum_Split
%The HighRe cases are anything equal to or above AOA_ReNum_Split
AOA_ReNum_Split = 1e6;
AOA_Low_LowRe = AOA_Low;
AOA_High_LowRe = AOA_High;
AOA_Low_HighRe = AOA_Low;
AOA_High_HighRe = AOA_High;

AOA_Low_Init = 0.0;
AOA_High = AOA_High;

AOA_delta = AOA_Delta;

MachInit = Mach;
ReInit = ReNumList(1);

%inout Low, Med High with 2 deltas
% or: [Low Del1 Med1 Del2 Med 2 Del3 High], and compute AOA vector which is
% then looked up later [-4 1 12 0.5 20 0.25 23]

%ReNumList = [.8E6; 1E6; 2E6; 3E6; 4E6; 6E6; 8E6; 10E6; 13E6];
%ReNumList = [1E6, ];
%ReNumList = [0.5E6; 0.6E6; 0.7e6; 0.8E6; 0.9E6 ];
%ReNumList = [10E6; 12e6; ];

       
 %% Start Java Robot and Click on RFOIL             
  rb = ROBOT_fnc;%sets 'rb' to stand for 'ROBOT_fnc' classdef
  rb.StartRobot()
                               
 disp('Start a new instance of RFOIL and Click on RFOIL window at root ("RFOIL   c>") level and DO NOT TOUCH ANYTHING UNTIL IT IS FINISHED');             
%give user a few seconds to click over to RFOIL
pause(rb.UserClickOverDelay)          

%make sure that we're at the RFOIL root
rb.PressEnter()
rb.PressEnter()
rb.PressEnter()
%% Load Airfoil Coordinates

if (0)
    %NACA 4 series airfoil
AirfoilNumber = '4412';
rb.LoadNACA4(AirfoilNumber);
else
    %Load Airfoil Coord. from file
    % AirfoilCoordFile = 'hand_18.dat';
    % AirfoilName = 'oso_18_hand1';
    rb.LoadAirfoil(AirfoilCoordFile,AirfoilName)
   
end

if (ChangePanel == 1)
    rb.PanelingChange(PanelNum);  %added 12-3-2024 dcm
end


%%  Make sure that ALFA_Temp_File does not exist 
%this is a temporary file used during each Alpha run to determine when each
%Alpha point is finished.  A new one is used for each Re# and deleted
%before the next Re#

%delete ALFA_Temp_File if it exists
if (exist(rb.ALFA_Temp_File));
delete(rb.ALFA_Temp_File);
end

%% Switch to analysis mode 'OPER'

rb.TypeStringPressEnter('OPER');
rb.TypeStringPressEnter('VISC');  % Switch to viscous mode, need to enter Reynolds Number next
   %set initial Reynolds Number
%ReInit = 1e6;
rb.TypeNumberPressEnter(ReInit);

%MachInit = 0.1;
rb.TypeStringPressEnter('Mach');
rb.TypeNumberPressEnter(MachInit);

    % START 6-20-13  Add initial run so that it doesn't crash when
    %trying to INIT an empty boudary layer...
    %calling INIT on the first run will cause aseq to crash.
    rb.TypeStringPressEnter('ALFA');
rb.TypeNumberPressEnter(AOA_Low_Init);
pause(rb.Pause_alfa)
 % Switch back to input promp
pause(rb.Pause_Keypress)
rb.SwitchToInputPrompt()%Switches from plot back to input prompt in Rfoil
pause(rb.Pause_Keypress)  
%END 6-20-13
% pause(1)
% disp('Hey You!  I''m Waiting! for you to click Matlab Commnd Window, Press Enter, and click back to RFoil prompt.')
% pause;  %added 6-13-13, waits for the user to click on the matlab prompt and press a key (ENTER works fine)
% pause(rb.UserClickOverDelay)%added 10-22-2024

%% Reynolds Number Loop

%input RENum list and alpha seq, could be multiple alpha range (fancy)
%create ReNum loop function in fnc, and create Alpha to call

for iRe = 1:length(ReNumList)
% iRe
ReNum = ReNumList(iRe)


%% Airfoil Loop Begin
numAOA = (AOA_High - AOA_Low_Init)/AOA_delta;
AOA_list = zeros(numAOA);
%Make AOA list
for iAOA = 1:(numAOA+1)  
    AOA_list(iAOA) = AOA_Low_Init + ((iAOA - 1) * AOA_delta);
end

%for iAOA = 1:(numAOA+1)
    
%AOA = AOA_Low + ((iAOA - 1) * AOA_delta)
    %AOA = 1.0;
    
    %Use this code to adjust the AOA limits for different Reynolds number:
    if (ReNum < AOA_ReNum_Split)
        AOA_Low = AOA_Low_LowRe;
        AOA_High = AOA_High_LowRe;
    else
        AOA_Low = AOA_Low_HighRe;
        AOA_High = AOA_High_HighRe;
    end


rb.TypeStringPressEnter('VPAR');
pause(rb.Pause_Keypress);%vpar needs an extra pause...
rb.TypeStringPressEnter('INIT');%calling INIT on the first run will cause aseq to crash.
rb.TypeStringPressEnter('RE');
rb.TypeNumberPressEnter(ReNum);

%dcm 10-23-2024 adding roughness/clean case
if RoughCleanFlag == 'C' %Clean case
    N_Factor_VPAR = N_factor_Clean;%n=9 is typically 'clean, 0.075 Turbulence)
    XTR_S = XTR_S_Clean;%S = suction surface (upper)
    XTR_P = XTR_P_Clean;%P = pressure Surface (lower)
elseif RoughCleanFlag == 'R' %Rough Case
    N_Factor_VPAR = N_factor_Rough;%n=9 is typically 'clean, 0.075 Turbulence)
    XTR_S = XTR_S_Rough;%S = suction surface (upper)
    XTR_P = XTR_P_Rough;%P = pressure Surface (lower)
else 
    disp('WARNING: RoughCleanFlag not equal to C or R in input table. Defaulting to Clean.')
    N_Factor_VPAR = N_factor_Clean;%n=9 is typically 'clean, 0.075 Turbulence)
    XTR_S = XTR_S_Clean;%XTR/c, S = suction surface (upper) transition location, 1.0 for clean. 0.0 for Leading Edge
    XTR_P = XTR_P_Clean;%XTR/c, P = pressure Surface (lower)
pause; 
end;
%Enter N and XTR values.
rb.TypeStringPressEnter('N');
rb.TypeNumberPressEnter(N_Factor_VPAR);
rb.TypeStringPressEnter('XTR');
rb.TypeNumberPressEnter(XTR_S);
rb.TypeNumberPressEnter(XTR_P);

%go up one level from VPAR menu
rb.PressEnterPause()

if DEBUG == 1;
    disp('Hey You!  I''m Waiting! for you to click Matlab Commnd Window, Press Enter, and click back to RFoil prompt.')
    pause;
    pause(rb.UserClickOverDelay)
end;

if (AOA_INIT_EACH_RENUM == 1)%This is only really needed if starting at really low AOA (say below -5 deg.) and if doing multiple ReNum with big steps
    % START 6-20-13  Add initial run to help get the run started.
    %TODO: Create ALFA function
    rb.TypeStringPressEnter('ALFA');
    rb.TypeNumberPressEnter(0);
    pause(rb.Pause_alfa)

    % Switch back to input promp
    pause(rb.Pause_Keypress)
    rb.SwitchToInputPrompt()%Switches from plot back to input prompt in Rfoil
    pause(rb.Pause_Keypress)

    rb.TypeStringPressEnter('ALFA');
    rb.TypeNumberPressEnter(-4);
    pause(rb.Pause_alfa)
    % Switch back to input promp
    pause(rb.Pause_Keypress)
    rb.SwitchToInputPrompt()%Switches from plot back to input prompt in Rfoil
    pause(rb.Pause_Keypress)
end
    rb.TypeStringPressEnter('ALFA');
rb.TypeNumberPressEnter(AOA_Low);
pause(rb.Pause_alfa)
 % Switch back to input promp
pause(rb.Pause_Keypress)
rb.SwitchToInputPrompt()%Switches from plot back to input prompt in Rfoil
pause(rb.Pause_Keypress)  
%END 6-20-13

%set write to temp file
rb.TypeStringPressEnter('PACC');
%This line outputs airfoilFileOut#.dat, where # is just an incremental
%numberof each Reynolds number file:
%ALFA_Out_File_i = strcat(AirfoilPolarOutFile,num2str(iRe),'.dat');
%This line outputs airfoilFileOut_r0010e4.dat, where this example would represent 1e6 Reynolds number:
%ALFA_Out_File_i = strcat(AirfoilPolarOutFile,'_',num2str(ReNum/(1.0e6)),'e6.dat');%Removed 5a 12-3-2024 to directly output desired filename and use tempfile here
ALFA_Out_File_i = strcat(AirfoilPolarOutFile);%,'_',num2str(ReNum/(1.0e6)),'e6.dat'); %Added 5a 12-3-2024 to directly output desired filename and use tempfile here
rb.TypeStringPressEnter(ALFA_Out_File_i);

%the following is in case a dump file exists, if it doens't than no harm
%done, and it will return to 'OPER'
%START 6-20-13 having issue with a crash in RFoil when the ASEQ analysis
%starts.  Try to fix by removing the following, since the same case works
%when I run it manually
rb.GetPastDumpFile()
%This wasn't the problem, INIT was the problem; although this really isn't
%needed anymore...
%rb.PressEnterPause();%in case there is a dump file there;%added 6-20-13, assumes dump file does not already exist
%END 6-20-13

%% Alpha analysis
% rb.TypeStringPressEnter('ALFA');
% rb.TypeNumberPressEnter(AOA);

rb.TypeStringPressEnter('ASEQ');
rb.TypeNumberPressEnter(AOA_Low);
rb.TypeNumberPressEnter(AOA_High);
rb.TypeNumberPressEnter(AOA_delta);

%pause(1.3)%replace with function to check if alfa run finished

%% Check if RFOIL airfoil polar File written

% START added 6-13-13
% RFoil proved to be unstable and crash for unknown reasons when using this
% scripting tool
%One guess is that the problem is the CheckPolarWritten function, which
%tries to read the output file as RFoil is writing to it, but this has not
%been proven.
%An attempted fix is to remove this function and replace it with a user
%return-press in matlab (pause), followed by a pause to allow the operator
%to click back to the RFoil window

%START 7-19-13 
if USER_CLICKOVER_INPUT == 1
%rb.CheckPolarWritten_ASEQ(ALFA_Out_File_i, AOA_list); %removed 6-13-13
disp('Hey You!  I''m Waiting! for you to click Matlab Commnd Window, Press Enter, and click back to RFoil prompt.')
pause;  %added 6-13-13, waits for the user to click on the matlab prompt and press a key (ENTER works fine)
pause(rb.UserClickOverDelay)%added 6-13-13
% END added 6-13-13
else
   disp(sprintf('Delay of %f seconds.', AOASEQ_DELAY)); %amount of delay to allow AOA loop to run, most runs are about 20 seconds, this should be very conservative
pause(AOASEQ_DELAY);
end
%END 7-19-13 

 %% Switch back to input promp
pause(rb.Pause_Keypress)
rb.SwitchToInputPrompt()%Switches from plot back to input prompt in Rfoil

%%

pause(rb.Pause_Keypress)         

% 	Pfil
% AOA_out_file.dat
% rb.TypeStringPressEnter('PFIL');
% rb.TypeStringPressEnter(strcat(AirfoilPolarOutFile,num2str(iRe),'.dat'));
% rb.GetPastDumpFile()
% % Padd
% rb.TypeStringPressEnter('PADD');
% % Pacc, disable auto point accumulation
rb.TypeStringPressEnter('PACC');

%% END AOA Loop
%end
%% END Reynolds Number Loop


%%  Make sure that ALFA_Temp_File does not exist 
%this is a temporary file used during each Alpha run to determine when each
%Alpha point is finished.  A new one is used for each Re# and deleted
%before the next Re#
%delete Aoa_temp1.dat, only delete after you're sure that the output file
%has been written!

%START 6-13-13, added check to see if the file exists before deleting it
  %note, this may be covering up other issues!  dcm...
%delete ALFA_Temp_File if it exists
if (exist(rb.ALFA_Temp_File));
delete(rb.ALFA_Temp_File);
end
%END 6-13-13, added check to see if the file exists before deleting it

end

%%  Output runtime

% TimeEndDays = now - TimeStartDays;
% RunTimeSeconds = TimeEndDays * (24.0 * 60.0);%convert serial time in days to minutes
% disp(RunTimeSeconds)

%%  End by going back to runnable state
%Press enter to root RFOIL level
pause(1.3)%wait for things to wrap up
for i = 1:4
    rb.PressEnter()
end
rb.PressEnterPause()
%return to inviscid so read for next call
rb.TypeStringPressEnter('OPER');
rb.TypeStringPressEnter('VISC');  % Switch to inviscous mode, need to enter Reynolds Number next
rb.PressEnter()

end