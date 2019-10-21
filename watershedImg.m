%Parameters for pre-processing and watershed.
addpath('Autosegmentation/')


frameFiles = dir('data\output\Position_*');
segmentation(frameFiles)

%% Tracking by taking into account the centroid and the overlapping between the previous

%Search all the images at 'data/output/segmentedcells/Image*'
watershededFiles = dir ('data/output/segmentedcells/Image #*');
figure;
for timepoint= 1:length(watershededFiles)-1
    timepoint
    actualTimePointImg = imread(strcat(watershededFiles(timepoint).folder, '/', watershededFiles(timepoint).name));
    centroidInfo= regionprops (actualTimePointImg, 'centroid');
    nextTimePointImg = imread (strcat(watershededFiles(timepoint+1).folder, '/', watershededFiles(timepoint+1).name));
    nextTimePointImg_labelled = bwlabel(nextTimePointImg);
    centroidIMG= cat (1, centroidInfo.Centroid);
    y = centroidIMG (:,2);
    x = centroidIMG(:,1);
    %     imshow(nextTimePointImg)
    %     hold on;
    pixelvalue = [];
     NewImg= zeros(size(nextTimePointImg_labelled),'like',nextTimePointImg_labelled);
    for numCell = 1:length(centroidInfo)
        %         plot(x(numCell), y(numCell), 'rx');
        pixelvalue(numCell, 1) = numCell;
        pixelvalue(numCell, 2) = nextTimePointImg_labelled(round(y(numCell)), round(x(numCell)));
        if pixelvalue (numCell, 2) ~=0
           
            
            NewImg(nextTimePointImg_labelled == pixelvalue(numCell, 1)) = pixelvalue(numCell, 2);
           
        end 
         
    end
    baseFileName = sprintf('Newimage #%03d.png', timepoint);
            Sust_Files= fullfile('data\output\NewImages',baseFileName);
            imwrite (NewImg, Sust_Files)
end


%     hold off

