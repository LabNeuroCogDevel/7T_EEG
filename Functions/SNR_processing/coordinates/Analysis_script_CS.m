
% Scrivener 23/02/21

% Script to open electrode coordinate files, convert to MNI space, and
%   save into a group results structure called 'Electrode_pleacement_results'

% The data is stored in the following folder structure:
%    ...manuscript_data\Participant_number
%       Participant_number = S2\S3\S4 etc.

% Within each participant folder is a text file with the electrode coordinates
%   eg: S2_electrode.txt
% As well as a text file with the snapped coordinates
%   eg: S2_electrodes_snapped.txt
% And finally the segmentation workspace with the affine transformation from SPM
%   eg: sEBRCB1_2-0006-00001-000208-01_seg8.mat

% Requirements:
%   The toolbox SPM
%   The function 'electrode2MNI.m' 
%   The function 'mni2orFROMxyz.m' 
%   The toolbox 'LoadNIIToolbox' should be in the path 
%   An empty nii file called 'EmptyMNI.nii' 




%_____________________________________________________________________________________________________

% LOOP 1

% Adds and saves the path to the main folder and subfolders
Directory_path = '...manuscript_data';
addpath(genpath(Directory_path));
Data_path = '...manuscript_data/electrode_files';


% File containing the cap size information
CapSize_FileName = ('CapSizes.txt');
CapSizes = readtable(CapSize_FileName);

% If the results structure doesn't already exist, create it 
if ~exist(fullfile(Data_path, 'Electrode_placement_results'))
   Electrode_placement_results = struct();
   Electrode_placement_results.Participant_number = 0;
   Electrode_placement_results.data = 0;
   cd(Directory_path)
   save('Electrode_placement_results','Electrode_placement_results')
end

% Declare the file suffixes
Electrode_suffix = '_electrodes.txt';
Snapped_suffix = '_electrodes_snapped.txt';


% Navigate to the data folder if not already there
cd (Data_path)

% Find the names of the participant folders
Ptp_folders = dir(Data_path);
% Delete the empty and random fields 
Ptp_folders = Ptp_folders(~ismember({Ptp_folders.name},{'.','..'}));
Ptp_folders = Ptp_folders(~ismember({Ptp_folders.name},{'.DS_Store','._.DS_Store'}));

% Reset the participant counter
participant_c = 1;    
    
% Loop over participants
for participant_c = 1:size(Ptp_folders,1)
    Folder_name = Ptp_folders(participant_c).name;
    Participant_number = Folder_name(2:end);

    % Navigate to the participants's folder
    Ptp_path = fullfile([Data_path '/' Folder_name]);
    cd (Ptp_path)

    % Load the electrode coordinate file
    Electrode_file_name = [num2str(Folder_name) Electrode_suffix];
    Electrode_data = readtable(Electrode_file_name);    

    % Load the snapped coordinate file
    Snapped_file_name = [num2str(Folder_name) Snapped_suffix];
    Snapped_data = readtable(Snapped_file_name);

     % Sort the data by electrode name
    Electrode_data = sortrows(Electrode_data,1);
    Snapped_data = sortrows(Snapped_data,1);

    % Load the affine matrix for the participant
    %   (Affine matrix is strored in Affine_workspace.Affine)
    Affine_file_name = ['sEBRCB1_' num2str(Participant_number) '-*.mat'];
    Affine_file = dir(Affine_file_name);
    Affine_workspace = load(Affine_file.name);

    % Copy the data matricies before editing
    Electrode_data_MNI = Electrode_data;
    Snapped_data_MNI = Snapped_data;

    % Pick out the X,Y,Z coordinates and open as an array
    Electrode_array = [Electrode_data.Loc_X, Electrode_data.Loc_Y, Electrode_data.Loc_Z];
    Snapped_array = [Snapped_data.Loc_X, Snapped_data.Loc_Y, Snapped_data.Loc_Z];

    % Create empty matricies for the results
    Electrode_MNI = nan(size(Electrode_array));
    Snapped_MNI = nan(size(Snapped_array));

    % Loop over electrodes
    % Convert the native coordinates to MNI space using function
    %   'electrode2MNI'
    for electrode_c = 1:size(Electrode_array,1)

        mni = electrode2MNI(Electrode_array(electrode_c,1), Electrode_array(electrode_c,2),...
              Electrode_array(electrode_c,3), Affine_workspace.Affine);
        Electrode_MNI(electrode_c,:) = mni(1:3);
        clear mni

        mni = electrode2MNI(Snapped_array(electrode_c,1), Snapped_array(electrode_c,2),...
              Snapped_array(electrode_c,3), Affine_workspace.Affine);
        Snapped_MNI(electrode_c,:) = mni(1:3);
        clear mni

    end

    clear Affine_workspace Affine_file Affine_file_name

    % Update the data files
    Electrode_data_MNI(:,2:4) = array2table(Electrode_MNI);
    Snapped_data_MNI(:,2:4) = array2table(Snapped_MNI);

    % Save the results into the results structure
    %load('Electrode_placement_results')
    Entry_number = size(Electrode_placement_results.data,2);
    Electrode_placement_results.data(Entry_number+1).Participant_number = Participant_number;
    Electrode_placement_results.data(Entry_number+1).Electrode = Electrode_data;
    Electrode_placement_results.data(Entry_number+1).Snapped = Snapped_data;
    Electrode_placement_results.data(Entry_number+1).Electrode_MNI = Electrode_data_MNI;
    Electrode_placement_results.data(Entry_number+1).Snapped_MNI = Snapped_data_MNI;

    % Save this in the main directory
    cd(Directory_path)
    save('Electrode_placement_results','Electrode_placement_results')

    % Go back to the rater folder to continue with the next subject
    cd(Data_path)
end

% Clear the empty field 
Electrode_placement_results.data(1) = [];


%_____________________________________________________________________________________________________

% Add cap sizes to Electrode_placement_results

% Declare the file name, load, and save
CapSize_FileName = ('CapSizes.txt');
Electrode_placement_results.CapSizes = readtable(CapSize_FileName);
cd(Directory_path)
save('Electrode_placement_results','Electrode_placement_results')


%_____________________________________________________________________________________________________

% LOOP 2

% Code to calculate the mean and SD of the electrode and brain locations across
%   participants (in MNI space), and add the results to
%   'Electrode_placement_results.mat'


% Navigate to the data folder

cd (Data_path)

% Find the names of the participant folders
Ptp_folders = dir(Data_path);
% Delete the empty and random fields 
Ptp_folders = Ptp_folders(~ismember({Ptp_folders.name},{'.','..'}));
Ptp_folders = Ptp_folders(~ismember({Ptp_folders.name},{'.DS_Store','._.DS_Store'}));

% Clear the variables
Electrode_MNI_XYZ = table;
Snapped_MNI_XYZ = table;

% Reset the participant number
participant_c = 1;

% Clear used variables
clear Brain_Average_X Brain_Average_Y Brain_Average_Z Elec_Average_X Elec_Average_Y Elec_Average_Z

for participant_c = 1:size(Electrode_placement_results.data,2)%size(Ptp_folders,1)
    Folder_name = Ptp_folders(participant_c).name;
    Participant_number = Folder_name(2:end);

    % Navigate to the participants's folder
    Ptp_path = fullfile([Data_path '/' Folder_name]);
    cd (Ptp_path)        

    % For each electrode, at X Y and Z separately, saves the values for
    %   all subjects into one table(electrode x subject)
    for electrode_c = 1:size(Electrode_placement_results.data(participant_c).Electrode,1)

       
        Elec_Average_X(electrode_c,participant_c) = Electrode_placement_results.data(participant_c).Electrode_MNI_XYZ(electrode_c,'Loc_X');
        Elec_Average_Y(electrode_c,participant_c) = Electrode_placement_results.data(participant_c).Electrode_MNI_XYZ(electrode_c,'Loc_Y');
        Elec_Average_Z(electrode_c,participant_c) = Electrode_placement_results.data(participant_c).Electrode_MNI_XYZ(electrode_c,'Loc_Z');

        Brain_Average_X(electrode_c,participant_c) = Electrode_placement_results.data(participant_c).Snapped_MNI_XYZ(electrode_c,'Loc_X');
        Brain_Average_Y(electrode_c,participant_c) = Electrode_placement_results.data(participant_c).Snapped_MNI_XYZ(electrode_c,'Loc_Y');
        Brain_Average_Z(electrode_c,participant_c) = Electrode_placement_results.data(participant_c).Snapped_MNI_XYZ(electrode_c,'Loc_Z');

    end

% Use the participant number to name the table columns        
Elec_Average_X.Properties.VariableNames(participant_c) = {Folder_name};
Elec_Average_Y.Properties.VariableNames(participant_c) = {Folder_name};
Elec_Average_Z.Properties.VariableNames(participant_c) = {Folder_name};

Brain_Average_X.Properties.VariableNames(participant_c) = {Folder_name};
Brain_Average_Y.Properties.VariableNames(participant_c) = {Folder_name};
Brain_Average_Z.Properties.VariableNames(participant_c) = {Folder_name};

end

% Add additional collumns for the mean and SD of locations
Elec_Average_X.mean = mean(table2array((Elec_Average_X)),2); Elec_Average_X.SD = std(table2array(Elec_Average_X),0,2);
Elec_Average_Y.mean = mean(table2array((Elec_Average_Y)),2); Elec_Average_Y.SD = std(table2array(Elec_Average_Y),0,2);
Elec_Average_Z.mean = mean(table2array((Elec_Average_Z)),2); Elec_Average_Z.SD = std(table2array(Elec_Average_Z),0,2);

Brain_Average_X.mean = mean(table2array((Brain_Average_X)),2); Brain_Average_X.SD = std(table2array(Brain_Average_X),0,2);
Brain_Average_Y.mean = mean(table2array((Brain_Average_Y)),2); Brain_Average_Y.SD = std(table2array(Brain_Average_Y),0,2);
Brain_Average_Z.mean = mean(table2array((Brain_Average_Z)),2); Brain_Average_Z.SD = std(table2array(Brain_Average_Z),0,2);

% Save them into the results structure
Electrode_placement_results.Averages.data.Elec_Average_X = Elec_Average_X;
Electrode_placement_results.Averages.data.Elec_Average_Y = Elec_Average_Y;
Electrode_placement_results.Averages.data.Elec_Average_Z = Elec_Average_Z;

Electrode_placement_results.Averages.data.Brain_Average_X = Brain_Average_X;
Electrode_placement_results.Averages.data.Brain_Average_Y = Brain_Average_Y;
Electrode_placement_results.Averages.data.Brain_Average_Z = Brain_Average_Z;



% Save the results structure in the main directory
cd(Directory_path)
save('Electrode_placement_results','Electrode_placement_results')

%_____________________________________________________________________________________________________
% LOOP 3

% Code to create masks for electrode and brain co-ordinates


% Add path to the load nii toolbox


% Load the template MNI nii file from SPM
EmptyMNI = load_nii('EmptyMNI.nii');

% Find the image dimentions
Image_size = size(EmptyMNI.img);

% Find the voxel size
Voxel_size = max(EmptyMNI.hdr.dime.pixdim(2:4));

% Create spheres for plotting based on euclidian distance approach
sphere_radius = 3; %in mm
r = sphere_radius -1;
t = randn((r*2+1),(r*2+1),(r*2+1));
dist = nan(1,numel(t));
center = ceil(size(t)/2) ;
for ii = 1:numel(t)
    [dd(1,1),dd(1,2),dd(1,3)]=ind2sub(size(t),ii);
    dd(2,:)=center;
    dist(ii)=pdist(dd)*mean(Voxel_size);
end
distmatrix = reshape(dist,size(t));
spheremask=zeros(size(t));
spheremask(distmatrix<=sphere_radius)=1;


% Converts the MNI co-ordinates to XYZ so that they can be passed to nii.img
for participant_c = 1:size(Electrode_placement_results.data,2)
    
    % Identify the participant number of the first entry 
    Participant_number = Electrode_placement_results.data(participant_c).Participant_number;
    Folder_name = ['S' Participant_number];
        
    % Navigate to the participants's folder
    Ptp_path = fullfile([Data_path '/' Folder_name]);
    cd (Ptp_path)
    
    % Clear used variables
    clear Brain_Average_X Brain_Average_Y Brain_Average_Z Elec_Average_X Elec_Average_Y Elec_Average_Z
    clear Electrode_MNI_XYZ Snapped_MNI_XYZ
    clear brain_nii brain_nii_conv elec_nii elec_nii_conv
    Electrode_MNI_XYZ = table;
    Snapped_MNI_XYZ = table;     
    
    % Create nii files based on the MNI default from SPM
    brain_nii = load_nii('EmptyMNI.nii');
    elec_nii = load_nii('EmptyMNI.nii');
    brain_nii_conv = load_nii('EmptyMNI.nii');
    elec_nii_conv = load_nii('EmptyMNI.nii');
    
    brain_nii_SD = load_nii('EmptyMNI.nii');
    brain_nii_SD_conv = load_nii('EmptyMNI.nii');
    
       
    
    for electrode_c = 1:size(Electrode_placement_results.data(participant_c).Electrode,1)
        
        % Convert Electrode_MNI to XYZ using function 'mni2orFROMxyz'
        [Electrode_MNI_XYZ.Loc_X(electrode_c,1), Electrode_MNI_XYZ.Loc_Y(electrode_c,1), Electrode_MNI_XYZ.Loc_Z(electrode_c,1)] = ...
            mni2orFROMxyz(Electrode_placement_results.data(participant_c).Electrode_MNI.Loc_X(electrode_c),...
            Electrode_placement_results.data(participant_c).Electrode_MNI.Loc_Y(electrode_c),...
            Electrode_placement_results.data(participant_c).Electrode_MNI.Loc_Z(electrode_c),Voxel_size,'mni');
        
        % Convert Snapped_MNI to XYZ using function 'mni2orFROMxyz'
        [Snapped_MNI_XYZ.Loc_X(electrode_c,1), Snapped_MNI_XYZ.Loc_Y(electrode_c,1), Snapped_MNI_XYZ.Loc_Z(electrode_c,1)] = ...
            mni2orFROMxyz(Electrode_placement_results.data(participant_c).Snapped_MNI.Loc_X(electrode_c),...
            Electrode_placement_results.data(participant_c).Snapped_MNI.Loc_Y(electrode_c),...
            Electrode_placement_results.data(participant_c).Snapped_MNI.Loc_Z(electrode_c),Voxel_size,'mni');
        
    end
    
    % Save the XYZ values into the results structure
    Electrode_placement_results.data(participant_c).Electrode_MNI_XYZ = Electrode_MNI_XYZ;
    Electrode_placement_results.data(participant_c).Snapped_MNI_XYZ = Snapped_MNI_XYZ;
        
    % Enter a value of 100 into the nii.img for each of the XYZ
    %   co-ordinates (the choice of 100 is arbitrary)  
    for electrode_c = 1:65
        brain_nii.img(round(Snapped_MNI_XYZ.Loc_X(electrode_c)),round(Snapped_MNI_XYZ.Loc_Y(electrode_c)),round(Snapped_MNI_XYZ.Loc_Z(electrode_c))) = 100;
        elec_nii.img(round(Electrode_MNI_XYZ.Loc_X(electrode_c)),round(Electrode_MNI_XYZ.Loc_Y(electrode_c)),round(Electrode_MNI_XYZ.Loc_Z(electrode_c))) = 100;
    end
    
    % Enter the SD into the nii.img for each of the XYZ
    %   co-ordinates     
    % First calculate mean SD for each electrode
    Electrode_placement_results.Averages.data.Brain_SD_average = ...
        mean([Electrode_placement_results.Averages.data.Brain_Average_X.SD,...
              Electrode_placement_results.Averages.data.Brain_Average_Y.SD,...
              Electrode_placement_results.Averages.data.Brain_Average_Z.SD],2);   
    for electrode_c = 1:65       
        brain_nii_SD.img(round(Snapped_MNI_XYZ.Loc_X(electrode_c)),round(Snapped_MNI_XYZ.Loc_Y(electrode_c)),round(Snapped_MNI_XYZ.Loc_Z(electrode_c))) = Electrode_placement_results.Averages.data.Brain_SD_average(electrode_c);
    end

        
    % Convolve the nii.img with spheres
    brain_nii_conv.img = convn(brain_nii.img,spheremask,'same');
    elec_nii_conv.img = convn(elec_nii.img,spheremask,'same');
    brain_nii_SD_conv.img = convn(brain_nii_SD.img,spheremask,'same');
   
    
    % Save the nii masks for the participant        
    save_nii(brain_nii, [Folder_name '_Brain_mean.nii']);
    save_nii(brain_nii_conv, [Folder_name '_Brain_mean_conv.nii']);
    save_nii(brain_nii_SD_conv, 'Brain_SD_conv3.nii');
 
    save_nii(elec_nii, [Folder_name '_Elec_mean.nii']);
    save_nii(elec_nii_conv, [Folder_name '_Elec_mean_conv.nii']);
    
    % Go back to the data folder to continue with the next subject
    cd(Data_path)
end


% Save the results structure in the main directory
cd(Directory_path)
save('Electrode_placement_results','Electrode_placement_results')

%____________________________________________________________________

% Code to save main results to an excel document
results_table = table();
results_table.electrode_names = Electrode_placement_results.data(1).Electrode.x_TargetName;
results_table.elec_average_x  = Electrode_placement_results.Averages.data.Elec_Average_X.mean;
results_table.elec_SD_x       = Electrode_placement_results.Averages.data.Elec_Average_X.SD;
results_table.elec_average_y  = Electrode_placement_results.Averages.data.Elec_Average_Y.mean;
results_table.elec_SD_y       = Electrode_placement_results.Averages.data.Elec_Average_Y.SD;
results_table.elec_average_z  = Electrode_placement_results.Averages.data.Elec_Average_Z.mean;
results_table.elec_SD_z       = Electrode_placement_results.Averages.data.Elec_Average_Z.SD;

results_table.brain_average_x  = Electrode_placement_results.Averages.data.Brain_Average_X.mean;
results_table.brain_SD_x       = Electrode_placement_results.Averages.data.Brain_Average_X.SD;
results_table.brain_average_y  = Electrode_placement_results.Averages.data.Brain_Average_Y.mean;
results_table.brain_SD_y       = Electrode_placement_results.Averages.data.Brain_Average_Y.SD;
results_table.brain_average_z  = Electrode_placement_results.Averages.data.Brain_Average_Z.mean;
results_table.brain_SD_z       = Electrode_placement_results.Averages.data.Brain_Average_Z.SD;

% Save the results structure in the main directory
cd(Directory_path)
writetable(results_table,'results_table.csv')



%_____________________________________________________________________________________________________






