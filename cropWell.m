function [] = cropWell(experimentDir)
%CROP_WELL Summary of this function goes here
%   Detailed explanation goes here
%% Process images to get only the well and crop it from all the images
    frameFiles = dir(fullfile(experimentDir, 'Position_*'));
    %gfpFiles= dir('data/gfp/gfp_*');
    outputDir = strrep(experimentDir, 'RawData', 'Output/Cropwell');
    mkdir(outputDir);
    %mkdir(fullfile(gfpFiles(1).folder, 'output'));
    allBBox = [];
    for timePoint = 1:length(frameFiles) %% To crop manually change timepoint
         img_original = imread(fullfile(frameFiles(timePoint).folder, frameFiles(timePoint).name));
         %img_gfp= imread(fullfile(gfpFiles(timePoint).folder, gfpFiles(timePoint).name));
         img_bin = imbinarize(img_original);
         %circularPoints = imfill(imgradient(img_original{numFrame}/65536, 'central')>0, 'holes');
         %img_bin(circularPoints) = 0;
         img_bin_dilated = imdilate(img_bin, strel('disk', 2));
    %      img_bin_bridging = bwmorph(img_bin_dilated, 'bridge');
    %      img_bin_closed = imclose(img_bin_dilated, strel('disk', 1));
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
         %croppedgfp= imcrop(img_gfp, bbox.BoundingBox);
         imwrite(croppedImage, fullfile(outputDir, frameFiles(timePoint).name));
         %imwrite(croppedgfp, fullfile(gfpFiles(timePoint).folder, 'output', gfpFiles(timePoint).name));
    end
end

