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
previousTrackingCells = [];
for timepoint= 1:length(watershededFiles)-1
    timepoint
    %     actualTimePointImg = imread(strcat(watershededFiles(timepoint).folder, '/', watershededFiles(timepoint).name));
    %     centroidInfo= regionprops (actualTimePointImg, 'centroid');
    %     nextTimePointImg = imread (strcat(watershededFiles(timepoint+1).folder, '/', watershededFiles(timepoint+1).name));
    %     nextTimePointImg_labelled = bwlabel(nextTimePointImg);
    %     centroidIMG= cat (1, centroidInfo.Centroid);
    %     y = centroidIMG (:,2);
    %     x = centroidIMG(:,1);
    %     %     imshow(nextTimePointImg)
    %     %     hold on;
    %     pixelvalue = [];
    %      NewImg= zeros(size(nextTimePointImg_labelled),'like',nextTimePointImg_labelled);
    %     for numCell = 1:length(centroidInfo)
    %         %         plot(x(numCell), y(numCell), 'rx');
    %         pixelvalue(numCell, 1) = numCell;
    %         pixelvalue(numCell, 2) = nextTimePointImg_labelled(round(y(numCell)), round(x(numCell)));
    %         if pixelvalue (numCell, 2) ~=0
    %
    %
    %             NewImg(nextTimePointImg_labelled == pixelvalue(numCell, 1)) = pixelvalue(numCell, 2);
    %
    %         end
    %
    %     end
    %     baseFileName = sprintf('Newimage #%03d.png', timepoint);
    %             Sust_Files= fullfile('data\output\NewImages',baseFileName);
    %             imwrite (NewImg, Sust_Files)
    %             %find matches
    
    % matchedfeatures
    %     frameFiles = dir('data\output\segmentedcells\Image #*');
    I1 = imread(strcat(watershededFiles(timepoint).folder, '/', watershededFiles(timepoint).name));
    I2 = imread(strcat(watershededFiles(timepoint+1).folder, '/', watershededFiles(timepoint+1).name));
    % Find the corners.
    
    % Extract features from the images:
    % Possible methods:
    %     detectBRISKFeatures
    %     detectFASTFeatures
    %     detectKAZEFeatures
    %     detectHARRISFeatures
    %     detectMinEigenFeatures
    %     detectMSERFeatures
    %     detectSURFFeatures
    
    % This may be replaced by centroid function
    points1 = detectKAZEFeatures(I1);
    points2 = detectKAZEFeatures(I2);
    % Extract the neighborhood features.
    
    [features1,valid_points1] = extractFeatures(I1,points1, 'blocksize', 11);
    [features2,valid_points2] = extractFeatures(I2,points2, 'blocksize', 11);
    % Match the features.
    resizedImg1 = imresize(I1, max(size(I1), size(I2)))>0;
    resizedImg2 = imresize(I2, max(size(I1), size(I2)))>0;
    combinedImg = imfuse(resizedImg1, resizedImg2,'blend','Scaling','joint')>0;%common and differences in diff colours
    combinedImg2 = imfuse(resizedImg1, resizedImg2,'diff','Scaling','joint')>0;%just differences
    OverlappingImg= combinedImg - combinedImg2;
    labelledImg1 = bwlabel(resizedImg1);
    labelledImg2 = bwlabel(resizedImg2);
    
    %%
    mkdir(fullfile(frameFiles(1).folder, 'Tracking'));
    NewImg= zeros(size(labelledImg2),'like',labelledImg2);
    trackingCells = {};
    for numCell= 1:max (max (labelledImg1))
        numCell;
        
        trackingCells(numCell, 1) = {numCell};
        uniqueLabels = unique(labelledImg2(labelledImg1 == numCell));
        uniqueLabels(uniqueLabels == 0) = [];
        trackingCells(numCell, 2) = {uniqueLabels};
        
        if length(trackingCells{numCell, 2}) == 1
            if isempty(previousTrackingCells) == 0 && size(previousTrackingCells, 1) >= size(trackingCells, 1)
                %% CARE 
                trackingCells(numCell, 2) = {previousTrackingCells{trackingCells{numCell, 1}, 2}};
            end
        else
            timepoint
            numCell
            length(trackingCells{numCell, 2})
            dividingCells = trackingCells{numCell, 2};
            areas = regionprops(labelledImg2,'Area');
            [~, indexMother] = max([areas(trackingCells{numCell, 2}).Area]);
            [~, indexDaughter] = min([areas(trackingCells{numCell, 2}).Area]);
            trackingCells(numCell, 2) = {dividingCells(indexMother)};
            trackingCells(max (max (labelledImg1))+1, 2) = {dividingCells(indexDaughter)};
            trackingCells(max (max (labelledImg1))+1, 1) = {max(max(labelledImg1))+1};
        end
        NewImg(labelledImg2 == trackingCells{numCell, 1}) = trackingCells{numCell, 2}; 
    end
    previousTrackingCells = trackingCells;
    
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
