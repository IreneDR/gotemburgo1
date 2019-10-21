function [] = segmentation(frameFiles)
%SEGMENTATION Summary of this function goes here
%   Detailed explanation goes here
    mkdir(fullfile(frameFiles(1).folder, 'labelledCells'));
    mkdir(fullfile(frameFiles(1).folder, 'objectsToKeep'));
    mkdir(fullfile(frameFiles(1).folder, 'allObjects'));
    mkdir(fullfile(frameFiles(1).folder,'segmentedcells'));
    mkdir(fullfile(frameFiles(1).folder,'segmentedcells_coloured'));
    mkdir(fullfile(frameFiles(1).folder, 'NewImages'));

    for timepoint= 1:length(frameFiles)
        timepoint

        IB= imread(fullfile(frameFiles(timepoint).folder, frameFiles(timepoint).name));
        IB=double(IB)/65536;
        filledWell = imfill(imclose(imbinarize(IB, 'adaptive', 'ForegroundPolarity','dark'), strel('disk', 1)), 'holes');


        insideWell = bwmorph(filledWell, 'majority');
        %dilatedIB = imdilate (insideWell, strel ('disk',2));
        IB(imerode(bwareafilt(insideWell, 1), strel('disk', 20)) == 0) = 0;%%Da problemas. No rellena el pocillo en las 147-150

        %Gaussian filter
        IB_filt=imgaussfilt(IB,6);
        IB=IB-IB_filt;
        IB=imcomplement(IB);
        %Normalize 0-255
        %IB=((IB-min(IB(:)))./(max(IB(:))-min(IB(:)))).*255;
        %% Binarize
        %level = graythresh(IB)

        binarizedImg = imbinarize(IB, 'adaptive', 'sensitivity', 0.7)==0; %0.7 good results
        baseFileName = sprintf('Image #%d.png', timepoint);
        fullFileName = fullfile('data\output\objectsToKeep', baseFileName);
        imwrite (binarizedImg, fullFileName)

        %     binarizedImg = imbinarize(IB, 'adaptive', 'sensitivity', 0.67)==0;
        %     baseFileName = sprintf('Image #%d.png', timepoint);
        %     fullFileName = fullfile('data\output\test5', baseFileName);
        %     imwrite (binarizedImg, fullFileName)
        localThreshold = multithresh(IB);
        binarizedImg3 = imbinarize(IB, localThreshold * 1.02)==0;
        baseFileName3 = sprintf('Image #%d.png', timepoint);
        fullFileName3 = fullfile('data\output\allObjects', baseFileName3);
        imwrite (binarizedImg3, fullFileName3)

        labelledImg = bwlabel(binarizedImg3, 4);
        labelledImgOnlyCells = ismember(labelledImg, unique(labelledImg(binarizedImg)));
        baseFileName4 = sprintf('Image #%d.png', timepoint);
        fullFileName4 = fullfile('data\output\labelledCells', baseFileName4);
        imwrite (labelledImgOnlyCells, fullFileName4)

        %% Outline of the well
        outline = imdilate(bwareafilt(labelledImgOnlyCells, 1, 8), strel('disk', 2));
        labelledImgOnlyCells(outline) = 0;
        %% Fill cluster of cells
        %closeimg= imclose (labelledImgOnlyCells, strel('disk',1));
        closeimg= imfill (labelledImgOnlyCells, 'holes');
        %smallBW2= bwareafilt(closeimg == 0,[1 7], 4);
        %closeimg = closeimg + smallBW2;
        baseFileName = sprintf('ImageIBclosed #%d.png', timepoint);
        fullFileName = fullfile('data\output\segmentedcells', baseFileName);
        imwrite (closeimg, fullFileName)

        %closeimg= imfill (closeimg, 'holes');
        %     suma= double(labelledImgOnlyCells) + double(fillimg);
        %figure, imshow(suma+1, parula(3))
        %figure, imshow(suma == 1)

        %% Labelling cells and removing artifacts
        %     labelled= bwlabel (closeimg==1, 4);
        %     biggestObject = bwareafilt(closeimg>0, 1);
        %     closeimg(biggestObject) = 0;
        %figure, imshow(biggestObject)
        %% Removing small objects
        %smallBW2= bwareafilt(suma == 1,[1 7], 4);
        %figure, imshow (smallBW2)
        %figure, imshow (labeledsuma+1, parula (20))
        %     sumafinal= suma>0;
        %figure, imshow (sumafinal)
        D = bwdist(~closeimg);
        %figure, imshow(D)
        %imshow(D,[],'InitialMagnification','fit')
        D = -D;
        D(closeimg==0) = Inf;
        Iwatershed = watershed(D);
        Iwatershed(~closeimg) = 0;
        rgb = label2rgb(Iwatershed,'jet',[.5 .5 .5]);
        %figure
        %imshow(rgb,'InitialMagnification','fit')
        %title('Watershed transform of D')
        baseFileName = sprintf('Image #%03d.png', timepoint);
        fullFileName = fullfile('data\output\segmentedcells', baseFileName);
        imwrite (Iwatershed>0, fullFileName)
        fullFileName = fullfile('data\output\segmentedcells_coloured', baseFileName);
        imwrite (rgb, fullFileName)

        %% Postprocessing to fix incorrectly separated cells

    end
end

