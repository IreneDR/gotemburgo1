%% Process images to get only the well and crop it from all the images

frameFiles = dir('First_ExperimentImageSequence\Position_*');
img_original = {};
mkdir(fullfile(frameFiles(1).folder, 'output'));

allBBox = [];
for numFrame = 1:length(frameFiles)
     img_original{numFrame, 1} = imread(fullfile(frameFiles(numFrame).folder, frameFiles(numFrame).name));
     img_bin = imbinarize(img_original{numFrame});
     %circularPoints = imfill(imgradient(img_original{numFrame}/65536, 'central')>0, 'holes');
     %img_bin(circularPoints) = 0;
     img_bin_dilated = imdilate(img_bin, strel('disk', 2));
%      img_bin_bridging = bwmorph(img_bin_dilated, 'bridge');
%      img_bin_closed = imclose(img_bin_dilated, strel('disk', 1));
     img_onlyMiddleCircle = bwareafilt(bwmorph(img_bin_dilated, 'hbreak'), 1, 8);
     bbox = regionprops(img_onlyMiddleCircle, 'boundingbox');
     if isempty(allBBox) == 0
         difference = allBBox(numFrame-1, :) - bbox.BoundingBox(3:4);
         if sum(difference) > 10
             disp('Error');
             [~, bbox.BoundingBox] = imcrop(img_original{numFrame, 1});
         end
     end
     allBBox(numFrame, :) = bbox.BoundingBox(3:4);
     img_original{numFrame, 2} = imcrop(img_original{numFrame}, bbox.BoundingBox);
     imwrite(img_original{numFrame, 2}, fullfile(frameFiles(numFrame).folder, 'output', frameFiles(numFrame).name));
end

