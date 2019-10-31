%Parameters for pre-processing and watershed.
addpath('Autosegmentation/')
addpath('src/')

frameFiles = dir('data\output\Position_*');
% segmentation(frameFiles)

%% Fixing segmentation errors
% for numFrame = [87 142]
%     baseFileName = sprintf('Image #%03d.png', numFrame);
%     fullFileName = fullfile('data\output\segmentedcells', baseFileName);
%     imgToFix = imread(fullFileName);
%     figure, imshow(imgToFix)
% %     imgToFix(103, 94)
% 
%     imwrite(imgToFix, fullFileName);
% end

%% Tracking by taking into account the centroid and the overlapping between the previous

%Search all the images at 'data/output/segmentedcells/Image*'
watershededFiles = dir ('data/output/segmentedcells/Image #*');
timepoint = 1;
I1 = imread(strcat(watershededFiles(timepoint).folder, '/', watershededFiles(timepoint).name));
previousNewImg = bwlabel(I1);

%% INIT HeritageTable 
HeritageInfo = {};
for ID_Cell= 1:max (max (previousNewImg))
    HeritageInfo = vertcat(HeritageInfo, {ID_Cell, [], []});
end

for timepoint= 2:length(watershededFiles)-1
    timepoint
    I2 = imread(strcat(watershededFiles(timepoint).folder, '/', watershededFiles(timepoint).name));
    resizedImg1 = imresize(previousNewImg, max(size(previousNewImg), size(I2)), 'nearest');
    resizedImg2 = imresize(I2, max(size(previousNewImg), size(I2)))>0;
    combinedImg = imfuse(resizedImg1>0, resizedImg2,'blend','Scaling','joint')>0;%common and differences in diff colours
    combinedImg2 = imfuse(resizedImg1>0, resizedImg2,'diff','Scaling','joint')>0;%just differences
    OverlappingImg = combinedImg - combinedImg2;
    labelledImg1 = resizedImg1;
    labelledImg2 = bwlabel(resizedImg2);
    OverlappingImg = bwlabel(OverlappingImg);

    
    %%
    mkdir(fullfile(frameFiles(1).folder, 'Tracking'));
    NewImg= zeros(size(labelledImg2),'like',labelledImg2);
    foundCells = [];
    for numCell= 1:max (max (labelledImg1))
        numCell;
        
        uniqueLabels = unique(labelledImg2(labelledImg1 == numCell));
        uniqueLabels(uniqueLabels == 0) = [];
        
        foundCells = unique(vertcat(uniqueLabels, foundCells));
        
        if length(uniqueLabels) == 1
             NewImg(labelledImg2 == uniqueLabels) = numCell;
        elseif length(uniqueLabels) > 1
            dividingCells = uniqueLabels;
            areas = regionprops(labelledImg2,'Area');
            [~, indexMother] = max([areas(dividingCells).Area]);
            [~, indexDaughter] = min([areas(dividingCells).Area]);
            %% The mother cell get its real ID
            NewImg(labelledImg2 == dividingCells(indexMother)) = numCell;
            %% We create a new ID for the daughter cell
            newDividingCell = max (max (labelledImg2));
            %% Assign the new dividing cell to its new ID
            NewImg(labelledImg2 == dividingCells(indexDaughter)) = newDividingCell;
            
            %% Assign daughter to mother
            HeritageInfo(numCell, 2) = {vertcat(HeritageInfo{numCell, 2}, newDividingCell)};
            HeritageInfo(numCell, 3)= {vertcat(HeritageInfo{numCell, 3}, timepoint)};
            
            HeritageInfo = vertcat(HeritageInfo, {newDividingCell, [], []});
        else
            disp('ERRRRRRROR!');
        end
    end
    
    % If we have cells without a correspondance on the NewImg, we need to
    % put a label on the NewImg, with the anew divided cells
    if length (unique(NewImg)) < length (unique(labelledImg2))
        newDividingCell = max (max (labelledImg2));
        labelledImg2(ismember(labelledImg2, foundCells)) = 0;
        uniqueLabels = unique(labelledImg2);
        uniqueLabels(uniqueLabels == 0) = [];
        centroidInfoOfMothers = regionprops (NewImg, 'centroid');
        newlabels = (1: length(uniqueLabels)) + length(foundCells);
        for numNewCell = 1:length(uniqueLabels)
            NewImg(labelledImg2 == uniqueLabels(numNewCell)) = newlabels(numNewCell);
            
            %% Centroid region of daughter
            imageWithOnlyDaughter = labelledImg2 == uniqueLabels(numNewCell);
            DaughterCentroid= regionprops (imageWithOnlyDaughter, 'centroid');
            [ClosestToNewCell, closestCellIndex]= pdist2(vertcat(centroidInfoOfMothers.Centroid), vertcat(DaughterCentroid.Centroid),'euclidean', 'smallest', 1)
            %% Do more stuff about daughter/mother related to centroid or closeness
            
            numCell = closestCellIndex;
            HeritageInfo(numCell, 2) = {vertcat(HeritageInfo{numCell, 2}, newDividingCell)};
            HeritageInfo(numCell, 3)= {vertcat(HeritageInfo{numCell, 3}, timepoint)};
            
            HeritageInfo = vertcat(HeritageInfo, {newDividingCell, [], []});
        end
    end
    
    previousNewImg = NewImg;
    
    baseFileName = sprintf('trackedimg #%03d.png', timepoint);
    Sust_Files= fullfile('data\output\Tracking',baseFileName);
    
    figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off'), imshow(NewImg+1,  colorcube(20))
    hold on;
    centroid = regionprops(NewImg, 'Centroid');
    centroids =round(vertcat(centroid.Centroid));
    for cellid= 1:size(centroids, 1)
        text(centroids(cellid,1), centroids(cellid,2), {cellid},'color', 'white' )
    end
    print (Sust_Files,'-dpng', '-r300')
    
    close all
end
