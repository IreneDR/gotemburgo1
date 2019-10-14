%% Process images to get only the well and crop it from all the images

frameFiles = dir('data\input\Position_*');
gfpFiles= dir('data/gfp/gfp_*');
mkdir(fullfile(frameFiles(1).folder, 'output'));
mkdir(fullfile(gfpFiles(1).folder, 'output'));
allBBox = [];
for timePoint = 1:length(frameFiles)
     img_original = imread(fullfile(frameFiles(timePoint).folder, frameFiles(timePoint).name));
     img_gfp= imread(fullfile(gfpFiles(timePoint).folder, gfpFiles(timePoint).name));
     img_bin = imbinarize(img_original);
     %circularPoints = imfill(imgradient(img_original{numFrame}/65536, 'central')>0, 'holes');
     %img_bin(circularPoints) = 0;
     img_bin_dilated = imdilate(img_bin, strel('disk', 2));
%      img_bin_bridging = bwmorph(img_bin_dilated, 'bridge');
%      img_bin_closed = imclose(img_bin_dilated, strel('disk', 1));
     img_onlyMiddleCircle = bwareafilt(bwmorph(img_bin_dilated, 'hbreak'), 1, 8);
     bbox = regionprops(img_onlyMiddleCircle, 'boundingbox');
     if isempty(allBBox) == 0
         difference = allBBox(timePoint-1, :) - bbox.BoundingBox(3:4);
         if sum(difference) > 10
                
                disp(['Error ' num2str(timePoint)]);
             [~, bbox.BoundingBox] = imcrop(img_original);
         end
     end
     allBBox(timePoint, :) = bbox.BoundingBox(3:4);
     croppedImage = imcrop(img_original, bbox.BoundingBox);
     croppedgfp= imcrop(img_gfp, bbox.BoundingBox);
     imwrite(croppedImage, fullfile(frameFiles(timePoint).folder, 'output', frameFiles(timePoint).name));
       imwrite(croppedgfp, fullfile(gfpFiles(timePoint).folder, 'output', gfpFiles(timePoint).name));
end

