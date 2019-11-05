
experimentDir = 'D:\Irene\gotemburgo1\data\RawData\050619\serie1\';

% %% Cropwell
%cropWell(experimentDir)
% % Fix Crop manually
% frameFiles=  dir(fullfile(experimentDir, 'Position_*'));
%   for timePoint = [152]
%       
%          img_original = imread(fullfile(frameFiles(timePoint).folder, frameFiles(timePoint).name));
%         [~, bbox.BoundingBox] = imcrop(img_original);
%          croppedImage = imcrop(img_original, bbox.BoundingBox);
%          outputDir= strrep(experimentDir, 'RawData', 'Output/Cropwell');
%          imwrite(croppedImage, fullfile(outputDir, frameFiles(timePoint).name));
%   end
%% Segmentation
disp('------------------------ Segmentation ------------------------')
% sensitivity:
% 
%segmentation(experimentDir, 0.71, 1.02, 5);

% Fixing segmentation errors
% for numFrame = 147 %[87 142]
%     baseFileName = sprintf('Position_%03d.tif', numFrame);
%     outputDir = strrep(experimentDir, 'RawData', 'Output/SegmentedCells');
%     fullFileName = fullfile(outputDir, baseFileName);
%     imgToFix = imread(fullFileName);
%     figure, imshow(imgToFix)
%     %imgToFix(103, 94)
%     % Correct cell lost from frame to frame
%     goodFrame = numFrame - 1;% Frame that contains the cell I want.
%     baseFileName = sprintf('Position_%03d.tif', goodFrame);
%     outputDir = strrep(experimentDir, 'RawData', 'Output/SegmentedCells');
%     fullFileName = fullfile(outputDir, baseFileName);
%     imgGood = imread(fullFileName);
%     imgGood_labelled = bwlabel(imgGood);
%     figure, imshow(imgGood_labelled)
%     prompt = 'IDs missing: []'; %specify IDs missing separated with ,
%     idsToAdd = input(prompt);
%     imgGood_resized = imresize(imgGood_labelled, size(imgToFix), 'nearest');
%     idIWantToChange = -1;% Find the ID of the lost cell and add it to next image
%     imgToFix(ismember(imgGood_resized, idsToAdd)) = 1;
%     figure, imshow(imgToFix)
%     imwrite(imgToFix, fullFileName);
%     close all
% end
%% Tracking
disp('------------------------ Tracking ------------------------')
Tracking(experimentDir)


