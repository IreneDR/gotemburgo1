
experimentDir = 'D:\Irene\gotemburgo1\data\RawData\050619\serie2\';
gfpDir= 'D:\Irene\gotemburgo1\data\RawData\050619\serie2\gfp\';
% %% Cropwell
% cropWell(experimentDir, fullfile(experimentDir, 'gfp'));
% close all
% % Fix Crop manually
% frameFiles=  dir(fullfile(experimentDir, 'Position_*'));
% gfpFiles= dir(fullfile(gfpDir, 'Position_*'));

%   for timePoint = [20]
%       
%       img_original = imread(fullfile(frameFiles(timePoint).folder, frameFiles(timePoint).name));
%       img_gfp= imread(fullfile(gfpFiles(timePoint).folder, gfpFiles(timePoint).name));
%       [~, bbox.BoundingBox] = imcrop(img_original);
%    
%       croppedImage = imcrop(img_original, bbox.BoundingBox);
%       croppedgfp= imcrop(img_gfp, bbox.BoundingBox);
%                   
%          outputDir= strrep(experimentDir, 'RawData', 'Output/Cropwell');
%          outputDirGFP = strrep(gfpDir, 'RawData', 'Output/Cropwell');
%          imwrite(croppedImage, fullfile(outputDir, frameFiles(timePoint).name));
%          imwrite(croppedgfp, fullfile(outputDirGFP, gfpFiles(timePoint).name));
%   end
%% Segmentation
disp('------------------------ Segmentation ------------------------')
% sensitivity:
% 
segmentation(experimentDir, 0.71, 1.02, 5);

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
% disp('------------------------ Tracking ------------------------')
% Tracking(experimentDir)
%% GFP over BF
%Search all the images at 'experimentDir*'

% THRESHOLD = 30000;

inputDir= strrep(experimentDir, 'RawData', 'Output/Tracking')
labelledFiles = dir(fullfile(inputDir, 'Position_*.mat'));
outputDir= strrep(experimentDir, 'RawData', 'Output/GenerationTime')
mkdir (outputDir)
inputDirGFP = strrep(gfpDir, 'RawData', 'Output/Cropwell');
outputDirGFP = strrep(gfpDir, 'RawData', 'Output/Peaks');
mkdir (outputDirGFP)
gfpFiles= dir(fullfile(inputDirGFP, 'Position_*'));

 CellCycleInfo = {};
for timepoint= 1: length(labelledFiles)
    timepoint
    load(strcat(labelledFiles(timepoint).folder, '/', labelledFiles(timepoint).name));
    GFPimg= imread(strcat(gfpFiles(timepoint).folder, '/', gfpFiles(timepoint).name));
    GFPimg_resized2= imresize(GFPimg, size(NewImg), 'nearest');
    GFPimg_resized= imresize(GFPimg, size(NewImg), 'nearest');
    GFPimg_resized(imerode(NewImg==0, strel('disk', 10)))=0;
    meanPixelValue= mean(GFPimg_resized(GFPimg_resized>0));
    threshold = meanPixelValue*3;
    %figure, imshow(GFPimg_resized)
    GFPimg_resized(NewImg==0) = 0;
    peaks= GFPimg_resized > threshold;
    peaks_nosmallareas = bwareaopen(peaks, 1);
    answer = '';
    while isequal(answer, 'No') == 0
        h= figure ('visible', 'off');
        imshow(GFPimg_resized2);
        hold on
         set(h, 'units','normalized','outerposition',[0 0 1 1]);
        ax = get(h, 'Children');
        set(ax,'Units','normalized')
        set(ax,'Position',[0 0 1 1])
        %Look for the peaks
        centroids = regionprops(peaks_nosmallareas, 'Centroid');
        centroids = vertcat(centroids.Centroid);
        %[xs, ys] = find(peaks_nosmallareas);
        if size(centroids, 1)>0
            [ys] = centroids(:, 1);
            [xs] = centroids(:, 2);
            for numX = 1:length(xs)  
                % Go through all the xs and ys
                plot(ys(numX),  xs(numX), 'rx')
            end
        end
    
        peaks_labelled = bwlabel(peaks_nosmallareas);
        answer = 'No';
%         answer = questdlg('Do you want to add or remove any GFP point?', 'GFP', 'Add', 'Remove', 'No', 'No');
%     
%         if isequal(answer, 'No') == 0
%             points = impoint(ax);
%             newPoint = round(getPosition(points));
%             if isequal(answer, 'Add')
%                 peaks_nosmallareas(newPoint(2), newPoint(1)) = 1;
%             elseif isequal(answer, 'Remove')
%                 %% Remove all the peaks of the selected cell
%                 peaks_nosmallareas(peaks_labelled(newPoint(2), newPoint(1)) == peaks_labelled) = 0;
%             end
%         end
    end
    
    baseFileName = sprintf('Position_%03d.png', timepoint);
    peakFile= fullfile(outputDirGFP, baseFileName);
    saveas (h, peakFile)
    close all
    labelOfCells = unique (NewImg);
    labelOfCells(labelOfCells==0) = [];
    
    for cellID= labelOfCells'
        if length(CellCycleInfo) < cellID
            actualCellCycleInfo = 0;
        else
            actualCellCycleInfo = CellCycleInfo{cellID};
        end
        
        
        if any(peaks_nosmallareas(NewImg==cellID))
            G1_G2 = 1;
        else
            G1_G2 = 2;
        end
        
        if isequal(actualCellCycleInfo, 0) == 0
         	if (1-mod(length(actualCellCycleInfo), 2))+1 == G1_G2 %% if it is still on the same G, it remains on the same column
                G1_G2 = length(actualCellCycleInfo);
            else %% Otherwise, it will increase the column
                G1_G2 = length(actualCellCycleInfo)+1;
            end
        end
        
        if size(actualCellCycleInfo, 2) < G1_G2
            actualCellCycleInfo(G1_G2) = 0;
        end
        
        actualCellCycleInfo(G1_G2)= actualCellCycleInfo(G1_G2) + 5;
        CellCycleInfo{cellID} = actualCellCycleInfo;
        
    end
    
end
outputDir= strrep(experimentDir, 'RawData', 'Output/FinalExcel');
mkdir(outputDir)
FileName=strcat(outputDir, 'CellCycle.xls');
tablaCellCycle=cell2table (CellCycleInfo')
 writetable(tablaCellCycle, FileName)

