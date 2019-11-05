function [] = segmentation(experimentDir, sensitivity, additionalBinarizedValue, wellThickness)
%SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here

    inputDir = strrep(experimentDir, 'RawData', 'Output/Cropwell');
    frameFiles = dir(fullfile(inputDir, 'Position_*'));
    %gfpFiles= dir('data/gfp/gfp_*');
    outputDir = strrep(experimentDir, 'RawData', 'Output/SegmentedCells');
    mkdir(outputDir);
   

    for timepoint= 1:length(frameFiles)
        timepoint

        IB= imread(fullfile(frameFiles(timepoint).folder, frameFiles(timepoint).name));
        IB=double(IB)/65536;
        filledWell = imfill(imclose(imbinarize(IB, 'adaptive', 'ForegroundPolarity','dark'), strel('disk', 1)), 'holes');
        insideWell = bwmorph(filledWell, 'majority');
        %dilatedIB = imdilate (insideWell, strel ('disk',2));
        % wallThickness = 20;
        IB(imerode(bwareafilt(insideWell, 1), strel('disk', wellThickness)) == 0) = 0;
        %Gaussian filter
        IB_filt=imgaussfilt(IB, 6);
        IB=IB-IB_filt;
        IB=imcomplement(IB);
 
        %% Binarize
   
        %Initial value = 0.7
        binarizedImg = imbinarize(IB, 'adaptive', 'sensitivity', sensitivity)==0; %Sensitivity may vary with images. WE CAN CHANGE IT.
    
        localThreshold = multithresh(IB);
        % additionalBinarizedValue = 1.02
        binarizedImg3 = imbinarize(IB, localThreshold * additionalBinarizedValue)==0;
    

        labelledImg = bwlabel(binarizedImg3, 4);
        labelledImgOnlyCells = ismember(labelledImg, unique(labelledImg(binarizedImg)));
  

        %% Outline of the well
        outline = imdilate(bwareafilt(labelledImgOnlyCells, 1, 8), strel('disk', 2));
        labelledImgOnlyCells(outline) = 0;
        %% Fill cluster of cells
        %closeimg= imclose (labelledImgOnlyCells, strel('disk',1));
        closeimg= imfill (labelledImgOnlyCells, 'holes');
        D = bwdist(~closeimg); 
        D = -D;
        D(closeimg==0) = Inf;
        Iwatershed = watershed(D);
        Iwatershed(~closeimg) = 0;
        rgb = label2rgb(Iwatershed,'jet',[.5 .5 .5]);
        fullFileName = fullfile(outputDir, frameFiles(timepoint).name);
        imwrite (Iwatershed>0, fullFileName)
    

    end
end

