%Parameters for pre-processing and watershed.
addpath('Autosegmentation/')
addpath('src/')

frameFiles = dir('data\output\Position_*');
% segmentation(frameFiles)

%% Fixing segmentation errors
% for numFrame = [44 68 70 89 107 110 138 139]
%     baseFileName = sprintf('Image #%03d.png', numFrame);
%     fullFileName = fullfile('data\output\segmentedcells', baseFileName);
%     imgToFix = imread(fullFileName);
%     figure, imshow(imgToFix)
%
%
%
%     imwrite(imgToFix, fullFileName);
% end

%% Tracking by taking into account the centroid and the overlapping between the previous

%Search all the images at 'data/output/segmentedcells/Image*'
watershededFiles = dir ('data/output/segmentedcells/Image #*');
timepoint = 1;
I1 = imread(strcat(watershededFiles(timepoint).folder, '/', watershededFiles(timepoint).name));
previousNewImg = bwlabel(I1);
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
    if max(OverlappingImg(:)) > max(labelledImg2(:))
        disp('CARE: Possible overlapping issue!');
        OverlappingImg_unrealCell = bwareafilt(OverlappingImg>0, max(OverlappingImg(:)) - max(labelledImg2(:)), 'smallest');
        labelledImg2(OverlappingImg_unrealCell>0) = 0;
        labelledImg1(OverlappingImg_unrealCell>0) = 0;
    end
    
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
        else
            disp('ERRRRRRROR!');
        end
    end
    
    if length (unique(NewImg)) < length (unique(labelledImg2))
        newDividingCell = max (max (labelledImg2));
        labelledImg2(ismember(labelledImg2, foundCells)) = 0;
        uniqueLabels = unique(labelledImg2);
        uniqueLabels(uniqueLabels == 0) = [];
        uniqueLabels
        if length(uniqueLabels)>1
            disp('Error');
        end
        NewImg(labelledImg2 == uniqueLabels) = newDividingCell;
    end
    
    previousNewImg = NewImg;
    
    baseFileName = sprintf('trackedimg #%03d.png', timepoint);
    Sust_Files= fullfile('data\output\Tracking',baseFileName);
    imwrite (NewImg+1,  colorcube(20), Sust_Files)
    
    %     figure, imshow(I1)
    %     hold on;
    %     validsfeatures1= [];
    %     validPointsInsideCells1= [];
    %     for numPoint = 1:length(valid_points1)
    %         isInsideCell = I1 (round (valid_points1.Location(numPoint, 2)), round( valid_points1.Location(numPoint, 1)));
    %
    %         if isInsideCell==1
    %             validsfeatures1 = vertcat(validsfeatures1 ,features1(numPoint, :));
    %             validPointsInsideCells1 = vertcat(validPointsInsideCells1, numPoint);
    %             plot(valid_points1.Location(numPoint, 1), valid_points1.Location(numPoint, 2), 'rx')
    %         end
    %
    %     end
    %
    %     figure, imshow(I2)
    %     hold on;
    %     validsfeatures2= [];
    %     validPointsInsideCells2 = [];
    %     for numPoint = 1:length(valid_points2)
    %         isInsideCell = I2 (round (valid_points2.Location(numPoint, 2)), round( valid_points2.Location(numPoint, 1)));
    %
    %         if isInsideCell==1
    %             validsfeatures2 = vertcat(validsfeatures2 ,features2(numPoint, :));
    %             validPointsInsideCells2 = vertcat(validPointsInsideCells2, numPoint);
    % %             plot(valid_points2.Location(numPoint, 1), valid_points2.Location(numPoint, 2), 'rx')
    %         end
    %     end
    
    % MatchThreshold: useless
    % MaxRatio: increase value return more matches. To delete the furthest
    % matches, decrease the value.
    %     indexPairs = matchFeatures(validsfeatures1,validsfeatures2, 'matchThreshold', 100, 'MaxRatio',0.1, 'unique', true);
    
    % Retrieve the locations of the corresponding points for each image.
    
    %     matchedPoints1 = valid_points1(validPointsInsideCells1(indexPairs(:,1)),:);
    %     matchedPoints2 = valid_points2(validPointsInsideCells2(indexPairs(:,2)),:);
    %     % Visualize the corresponding points. You can see the effect of translation between the two images despite several erroneous matches.
    
    figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off' ); showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2);
    mkdir('data\output\features');
    baseFileName = sprintf('Newimage #%03d.png', timepoint);
    Sust_Files= fullfile('data\output\features',baseFileName);
    print (Sust_Files,'-dpng', '-r300')
    close all
end
