function [] = cropWell(experimentDir, gfpDir)
%CROP_WELL Summary of this function goes here
%   Detailed explanation goes here
%% Process images to get only the well and crop it from all the images
    frameFiles = dir(fullfile(experimentDir, 'Position_*'));
    gfpFiles= dir(fullfile(gfpDir, 'Position_*'));
    outputDir = strrep(experimentDir, 'RawData', 'Output/Cropwell');
    mkdir(outputDir);
    outputDirGFP = strrep(gfpDir, 'RawData', 'Output/Cropwell');
    mkdir(outputDirGFP);
    allBBox = [];
    for timePoint = 1:length(frameFiles) %% To crop manually change timepoint
         img_original = imread(fullfile(frameFiles(timePoint).folder, frameFiles(timePoint).name));
         
         img_gfp= imread(fullfile(gfpFiles(timePoint).folder, gfpFiles(timePoint).name));
         
         img_bin = imbinarize(img_original);
         img_bin_dilated = imdilate(img_bin, strel('disk', 2));
         img_onlyMiddleCircle = bwareafilt(bwmorph(img_bin_dilated, 'hbreak'), 1, 8);
         bbox = regionprops(img_onlyMiddleCircle, 'boundingbox');
%          [~, bbox.BoundingBox] = imcrop(img_original);
         if isempty(allBBox) == 0
             difference = allBBox(timePoint-1, :) - bbox.BoundingBox(3:4);
             if sum(difference) > 10
                disp(['Error ' num2str(timePoint)]);
                [~, bbox.BoundingBox] = imcrop(img_original, bbox.BoundingBox);
             end
         end
         allBBox(timePoint, :) = bbox.BoundingBox(3:4);
         croppedImage = imcrop(img_original, bbox.BoundingBox);
         croppedgfp= imcrop(img_gfp, bbox.BoundingBox);
         imwrite(croppedImage, fullfile(outputDir, frameFiles(timePoint).name));
         
         imwrite(croppedgfp, fullfile(outputDirGFP, gfpFiles(timePoint).name));
    end
end

