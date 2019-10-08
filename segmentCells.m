frameFiles = dir('data/output/Position_*');
mkdir ('data', 'result');
figure;
for numFrame = 1:length(frameFiles)
     img_original = imread(fullfile('data/output', frameFiles(numFrame).name));
     img_bin = not(imbinarize(img_original));
     img_bin_dilated = imdilate(img_bin, strel('disk', 2));
     imshow(img_bin_dilated)
end
   img_onlyMiddleCircle = bwareafilt(bwmorph(img_bin_dilated, 'hbreak'), 1, 8);
     bbox = regionprops(img_onlyMiddleCircle, 'boundingbox');
       