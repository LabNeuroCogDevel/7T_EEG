
addpath(('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/1. Recording & Stimulation of iEEG - pbelchps files/Image Reconstruction/customFcns'))
addpath(genpath('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/abelCode/ripple/Tools/neuroshare'))

datapath = ('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/sEEG_backup/sEEG_rawData/');
patientAnatPath = ('/Users/shanemckeon/Library/CloudStorage/OneDrive - University of Pittsburgh/1. Recording & Stimulation of iEEG - pbelchps files/Image Reconstruction/patients/');
savePath = ('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/sEEG_backup/sEEG_rawData/');
entropyPath = ('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/entropy/individualSubjectFiles/');


agefile = readtable('/Users/shanemckeon/Library/CloudStorage/OneDrive-UniversityofPittsburgh/PBE Lab/PBE_agefile.csv');

%load in all the files
setfiles0 = dir([datapath, 'P*/Rest/rest*ns2*']);
setfiles = {};

for epo = 1:length(setfiles0)
    setfiles{epo,1} = fullfile(datapath, setfiles0(epo).name); % cell array with EEG file names
    % setfiles = arrayfun(@(x) fullfile(folder, x.name), setfiles(~[setfiles.isdir]),folder, 'Uni',0); % cell array with EEG file names
end


for j = 1:length(setfiles0)

    idvalues(j,:) = (setfiles0(j).folder(106:112));

    if ~isfile([entropyPath idvalues(j,:) '_permEntropy.csv'])

        if isfile([savePath, idvalues(j,:), '/Rest/data.mat']) % is there a .mat file for the patient
            load([savePath, idvalues(j,:), '/Rest/data.mat']);
            oneMin = data(1:60000, :);
            dataDown = downsample(oneMin, 2);
            montage = readtable([savePath, idvalues(j,:), '/Rest/montage', idvalues(j,:), '.xlsx']);
           
            % Identify rows where the first column has NaN in the last two rows
            rowsToRemove = isnan(montage{end-1:end, 1});

            if any(rowsToRemove)
                % Remove the identified rows
                montage(end-1:end, :) = [];
                disp('Rows removed.');
            else
                disp('No rows to remove.');
            end


            for c = 1:size(data, 2)
                [Perm(c,:), Pnorm(c,:), cPE(c,:)] = PermEn(oneMin(:,c), 'm', 3, 'tau', 1);

                % Mobj = MSobject("SampEn");
                % [MSx, Ci(:,c)] = MSEn(dataDown(:,c), Mobj, 'Scales', 20, 'Methodx', 'coarse', 'RadNew', 0, 'Plotx', false);

            end

            permTable = array2table(Perm, 'VariableNames', {'perm1', 'perm2', 'perm3'});
            pnormTable = array2table(Pnorm, 'VariableNames', {'pnorm1', 'pnorm2', 'pnorm3'});
            cPETable = array2table(cPE, 'VariableNames', {'cPE1', 'cPE2'});

            subjectTable = horzcat(montage, permTable, pnormTable, cPETable);

            % Create a new column with subject ID repeated for every row
            subjectIDColumn = repmat(idvalues(j,:), size(subjectTable, 1),1);

            % Add the new column to the existing table
            subjectTable = [table(subjectIDColumn, 'VariableNames', {'Subject'}), subjectTable];

          
            subjectSavePath = [entropyPath idvalues(j,:) '_permEntropy.csv'];
            writetable(subjectTable, subjectSavePath)

        end

    end

    clear Perm
    clear permTable
    clear pnormTable
    clear cPETable
    clear Pnorm
    clear cPE
end