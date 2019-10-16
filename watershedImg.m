%
%%%%%%%%%%%%%%%%%%%%%   Automated_Seeding.m   %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective: segmenting the last image of an image time series
%--------------------------------------------------------------------------
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY.
%--------------------------------------------------------------------------
% For feedback and questions please contact me at:
% Nihal.Temamogullari@UTSouthwestern.edu
% nezgiwood@gmail.com
%--------------------------------------------------------------------------
%This script takes as input an image and segments cells.
%This script requires:
%                               pre_processing_and_watershed.m
%                               correct_segmentation.m
%                               segmentation_subroutine_for_seeding.m
%                               retain_largest_obj.m

clearvars; clc; close all;


%==========================================================================
%Specify Parameters

%Parameters for pre-processing and watershed.
addpath('Autosegmentation/')

frameFiles = dir('data\output\Position_*');
mkdir(fullfile(frameFiles(1).folder, 'test'));
mkdir(fullfile(frameFiles(1).folder, 'test3'));
mkdir(fullfile(frameFiles(1).folder, 'test4'));
mkdir(fullfile(frameFiles(1).folder, 'test5'));
mkdir(fullfile(frameFiles(1).folder,'data\output\segmentedcells'));

%Fair: 7 8 12 100
%Bad: 4 46 113

for timepoint= 1:length(frameFiles)
    timepoint
    size_padding=0; %specifies the size of padding.
    colony_boundary_bias=1; %Population boundary bias is added to the pixels at the boundary of the population. Population boundary is determined using the coarse foregound-background segmentation.
    size_strel_filters=3; %specifies the size of the structuring element used for average and standard deviation filters.
    cell_boundary_threshold=1; %specifies the threshold for designating a pixel as a boundary pixel.
    
    IB= imread(fullfile(frameFiles(timepoint).folder, frameFiles(timepoint).name));
    IB=double(IB)/65536;
    insideWell = bwmorph(imfill(imclose(imbinarize(IB), strel('disk', 1)), 'holes'), 'majority');
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
    fullFileName = fullfile('data\output\test', baseFileName);
    imwrite (binarizedImg, fullFileName)
            
%     binarizedImg = imbinarize(IB, 'adaptive', 'sensitivity', 0.67)==0;
%     baseFileName = sprintf('Image #%d.png', timepoint);
%     fullFileName = fullfile('data\output\test5', baseFileName);
%     imwrite (binarizedImg, fullFileName)
    localThreshold = multithresh(IB);
    binarizedImg3 = imbinarize(IB, localThreshold * 1.001)==0;
    baseFileName3 = sprintf('Image #%d.png', timepoint);
    fullFileName3 = fullfile('data\output\test3', baseFileName3);
    imwrite (binarizedImg3, fullFileName3)
    
    labelledImg = bwlabel(binarizedImg3, 4);
    labelledImgOnlyCells = ismember(labelledImg, unique(labelledImg(binarizedImg)));
    baseFileName4 = sprintf('Image #%d.png', timepoint);
    fullFileName4 = fullfile('data\output\test4', baseFileName4);
    imwrite (labelledImgOnlyCells, fullFileName4)
    
    %% Outline of the well
    outline = imdilate(bwareafilt(labelledImgOnlyCells, 1, 8), strel('disk', 2));
    labelledImgOnlyCells(outline) = 0;
    %% Fill cluster of cells

%      
    closeimg= imdilate (labelledImgOnlyCells, strel('disk',1));
       
     smallBW2= bwareafilt(closeimg == 0,[1 7], 4);
   closeimg = closeimg + smallBW2;
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
    baseFileName = sprintf('Image #%d.png', timepoint);
    fullFileName = fullfile('data\output\segmentedcells', baseFileName);
    imwrite (rgb, fullFileName)
    
    
    %After reading the bright-field image, IB:
    
    size_strel_bg_1=0; %specifies the size of the structuring element used for connecting the possible gaps on the cell bundaries during the coarse foreground-background segmentation step.
    size_strel_bg_2=20; %specifies the size of the structuring element used for closing at the step of coarse foreground-background segmentation.
    min_colony=100; %minimum number of pixels that will be considered as a population of cells.
    clean_BW=10; %connected components smaller than clean_BW are removed at the determination of cell boundary pixels.
    save_name='Seed1';
    
    %IB=double(IB); IB=((IB-min(IB(:)))./(max(IB(:))-min(IB(:)))).*255; % Phase image is converted to double and is scaled between 0 and 255.
    parametersh=[size_strel_filters cell_boundary_threshold size_padding min_colony size_strel_bg_1 size_strel_bg_2 colony_boundary_bias clean_BW]; %parameters required for the function pre_processing_and_watershed
    
    
    %Set parameters required for the function segmentation_subroutine_for_seeding.m. For a more detailed explanation of them see Doncic et al. 2013 PLOS ONE.
    threshold=0.2; %segmentation threshold
    phase_subtraction_factor  =3; %factor that determines how much of the phase image is subtracted from the segmentation.
    dist_modifier = 3; %factor modulating the penalty term that is applied to pixels that were not part of the cell in the previous segmentation.
    max_area_increase_per_tp=0.18; %maximal allowed percent area increase before switching to a higher segmentation threshold.
    nrep=5; %number of consecutive segmentations (repetitions). The segmentation converges to the correct cell area through these consecutive segmentations.
    %Set parameters required for automated correction & fine-tuning steps.
    min_cell_size=7; %minimum cell size
    max_cell_size=500; %maximum cell size
    cell_margin=5; %margin around the cell. This parameter is used to cut an area around the cell to create a subimage with the cell at the center.
    %==========================================================================
    %##Pre_processing and Watershed Step## (Figure 1B)
    %[Iwatershed] =pre_processing_and_watershed(IB,parametersh); %watershed result
    %==========================================================================
%     %##Automated Correction & Fine-Tuning Step##(Figure 1D)
%     no_obj=max(Iwatershed(:)); %number of objects
%     cells_to_segment=1:no_obj;
%     IPP=padarray(IB,[size_padding size_padding],'symmetric'); %padded phase image
%     IPP=uint8(IPP);
%     [x_size,y_size]=size(IPP); %size of the padded image.
%     Scores=zeros(x_size,y_size,no_obj,'uint8'); %This will be a reduction variable in the parfor loop. It has the 'cell scores'.
%     Scores2=zeros(x_size,y_size,no_obj,'uint8'); %This will be a reduction variable in the parfor loop. It has the cell locations. It also keeps track of divisions.
%     numcells=ones(1,no_obj); %This variable keeps track of divisions. Although we start with 'no_obj' many cells, some will get divided, some will be deleted and thus the final number of cells will be different.
%     
%     for index=1:length(cells_to_segment)
%         no_cell_of_interest=cells_to_segment(index); %the current cell
%         stats=regionprops(logical(Iwatershed==no_cell_of_interest),'Centroid','BoundingBox');
%         %##Cut an area around the cell (subimage)##
%         bbox=round(stats(1).BoundingBox);
%         %The following four lines make sure that the column/row numbers do not exceed the image size.
%         lower_x_limit=max(1,bbox(1)-cell_margin);
%         upper_x_limit=min(y_size,bbox(1)+bbox(3)+cell_margin);
%         lower_y_limit=max(1,bbox(2)-cell_margin);
%         upper_y_limit=min(x_size,bbox(2)+bbox(4)+cell_margin);
%         x_cn=lower_x_limit:upper_x_limit; %column numbers that make up the subimage. Note that when the image is visualized, column numbers correspond to x-axis coordinates, thus, the variable is named after x-coordinates.
%         y_cn=lower_y_limit:upper_y_limit; %row numbers that make up the subimage. Note that when the image is visualized, row numbers correspond to y-axis coordinates, thus, the variable is named after y-coordinates.
%         subimage_phase=double(IPP(y_cn,x_cn));
%         cell_of_interest_segmentation=(Iwatershed==no_cell_of_interest); %get the coarse segmentation that we got in the pre-processing and watershed step.
%         subimage_segmentation=double(cell_of_interest_segmentation(y_cn,x_cn));
%         
%         %##Watershed with multiple thresholds## This function checks whether the cell needs to be divided, i.e. checks for under-segmentation. Thus, it enables the correction of undersegmentation mistakes as illustrated in Figure 2B.
%         [new_segmentation]=check_divide(subimage_phase,subimage_segmentation);
%         
%         %============== find connected components ========================
%         %remove small components
%         CC=bwconncomp(new_segmentation,4); %connected components
%         numPixels = cellfun(@numel,CC.PixelIdxList); %number of pixels that make up each connected component
%         ccrem=find(numPixels<min_cell_size | numPixels>max_cell_size); % connected components to remove: Remove connected components that are less than the minimum cell size or greater than the maximum cell size.
%         if isempty(ccrem)==0 %if there is any connected component to remove
%             for j=ccrem
%                 new_segmentation(CC.PixelIdxList{j})=0;
%             end
%         end
%         %update numcells
%         numcells(1,index)=numcells(1,index)-1+CC.NumObjects-length(ccrem); %Update the number of cells in this subimage.
%         %##Divide?## (Blue box in Figure 1D)
%         if numcells(1,index)>1 %if there is a division
%             [x_sizeh,y_sizeh]=size(subimage_phase); %size of the subimage
%             Scoreh=zeros(x_sizeh,y_sizeh,numcells(1,index)); %Create a matrix that will hold the scores of each piece. 'h' is added to mark the variables that are only be used inside the loop/condition they are first defined.
%             counter=1;
%             for cellh=setdiff(1:CC.NumObjects,ccrem)%go over each piece (i.e. connected component) that is not part of ccrem
%                 subimage_segmentationh=zeros(size(new_segmentation));  %'h' is added to mark the variables that are only be used inside the loop/condition they are first defined.
%                 subimage_segmentationh(CC.PixelIdxList{cellh})=1; %this holds the segmentation of only the connected component of interest, i.e. cellh
%                 %get a bounding box around the connected component of interest
%                 stats=regionprops(logical(subimage_segmentationh),'Centroid','BoundingBox');
%                 bbox=round(stats(1).BoundingBox);
%                 %The following four lines make sure that the column/row numbers do not exceed the image size.
%                 lower_x_limit=max(1,bbox(1)-cell_margin);
%                 upper_x_limit=min(y_sizeh,bbox(1)+bbox(3)+cell_margin);
%                 lower_y_limit=max(1,bbox(2)-cell_margin);
%                 upper_y_limit=min(x_sizeh,bbox(2)+bbox(4)+cell_margin);
%                 x_cnh=lower_x_limit:upper_x_limit; %column numbers that make up the subimage
%                 y_cnh=lower_y_limit:upper_y_limit; %row numbers that make up the subimage
%                 subimage_phaseh=(subimage_phase(y_cnh,x_cnh));
%                 subimage_segmentationh=(subimage_segmentationh(y_cnh,x_cnh));
%                 %##Segmentation multiple times for each piece## This step refines the cell boundaries.
%                 [~,new_segmentationh]=segmentation_subroutine_for_seeding(subimage_phaseh,subimage_segmentationh,phase_subtraction_factor,dist_modifier,threshold,nrep,max_area_increase_per_tp);
%                 new_segmentationh=retain_largest_obj(new_segmentationh);
%                 Scoreh(y_cnh,x_cnh,counter)=new_segmentationh;
%                 counter=counter+1;
%             end
%             %##Check Overlaps between pieces##
%             cellsremoved=[]; %connected components/pieces to be removed
%             cellstorun=[];  %connected components/pieces that will be fed into segmentation_subroutine_for_seeding again
%             for kk=1:numcells(1,index) %First Cell
%                 for jj=kk+1:numcells(1,index) %Second Cell
%                     if ismember(kk,cellsremoved)==0 && ismember(jj,cellsremoved)==0 %Proceed if both cells are not removed before.
%                         first_cell_segmentation= Scoreh(:,:,kk)>0;
%                         second_cell_segmentation=Scoreh(:,:,jj)>0;
%                         Int=first_cell_segmentation.*second_cell_segmentation; %intersection of the two cells
%                         Intsize=sum(Int(:)); %size of the intersection
%                         first_cell_size=sum(first_cell_segmentation(:));
%                         second_cell_size=sum(second_cell_segmentation(:));
%                         if max(Intsize/first_cell_size, Intsize/second_cell_size)>0.25  %##Merge?## %If the intersection size is larger than 25% of the smaller cell size then we merge the two cells.
%                             numcells(1,index)=numcells(1,index)-1; %since there is a merging event, number of cells here goes down by 1.
%                             %merge cells by adding the smaller cell to the larger cell and removing the smaller cell.
%                             if first_cell_size>second_cell_size
%                                 cellsremoved=[cellsremoved jj];
%                                 cellstorun=[cellstorun kk];
%                                 %remove the second cell and add parts of the second cell that are not in the first cell to the first cell
%                                 Scoreh(:,:,kk)=Scoreh(:,:,kk)+(Scoreh(:,:,jj).*(~Int));
%                             else
%                                 %the reverse
%                                 cellsremoved=[cellsremoved kk];
%                                 cellstorun=[cellstorun jj];
%                                 Scoreh(:,:,jj)=Scoreh(:,:,jj)+(Scoreh(:,:,kk).*(~Int));
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             %Running merged cells again, i.e. feeding these cells into the segmentation_subroutine_for_seeding again to refine the cell boundaries.
%             if isempty(cellstorun)==0
%                 for jj=unique(cellstorun)
%                     subimage_segmentationh=(Scoreh(:,:,jj)>0);
%                     stats=regionprops(logical(subimage_segmentationh),'Centroid','BoundingBox');
%                     bbox=round(stats(1).BoundingBox);
%                     %The following four lines make sure that the column/row numbers do not exceed the image size.
%                     lower_x_limit=max(1,bbox(1)-cell_margin);
%                     upper_x_limit=min(y_sizeh,bbox(1)+bbox(3)+cell_margin);
%                     lower_y_limit=max(1,bbox(2)-cell_margin);
%                     upper_y_limit=min(x_sizeh,bbox(2)+bbox(4)+cell_margin);
%                     x_cnh=lower_x_limit:upper_x_limit; %column numbers that make up the subimage
%                     y_cnh=lower_y_limit:upper_y_limit; %row numbers that make up the subimage
%                     subimage_phaseh=(subimage_phase(y_cnh,x_cnh));
%                     subimage_segmentationh=(subimage_segmentationh(y_cnh,x_cnh));
%                     %##Segmentation multiple times##
%                     [~,new_scoreh]=segmentation_subroutine_for_seeding(subimage_phaseh,subimage_segmentationh,phase_subtraction_factor,dist_modifier,threshold,nrep,max_area_increase_per_tp);
%                     new_scoreh=retain_largest_obj(new_scoreh);
%                     Scoreh(y_cnh,x_cnh,jj)=new_scoreh; % add new cell score to Scoreh
%                 end
%             end
%             
%             %Remove the scores of removed cells
%             for kk= sort(cellsremoved,'descend')
%                 Scoreh(:,:,kk)=[];
%             end
%             
%             %##Distribute overlaps using scores##
%             [~,Cell_Numbers] = max(Scoreh,[],3); %compare the score of every cell at every pixel and return the cell number with the highest score at each pixel.
%             %Note that the above line assigns 1 to the pixels where the first cell has the maximum score, but also to the background. To correct that we use the next two lines, which assign zero to background pixels.
%             BG=(Scoreh(:,:,1)==0).*(~(Cell_Numbers>1)); %background
%             Cell_Numbers(logical(BG))=0; %label matrix, where the location of each cell is marked by its number.
%             
%             for ii=1:numcells(1,index)
%                 cell_score=Scoreh(:,:,ii); %score of the current cell
%                 cell_score(~(Cell_Numbers==ii))=0; %remove the pixels from that cell score, where other cells have higher scores
%                 Scoreh(:,:,ii)=cell_score;
%             end
%             
%             new_segmentation=zeros(size(subimage_phase));
%             new_score=zeros(size(subimage_phase));
%             
%             %Note that the matrices Scoreh and Scoreh2 are created for this particular if condition ( if numcells(1,index)>1 %if there is a division). Now we put the information in them into the matries new_score and new_segmentation.
%             for ii=1:numcells(1,index)
%                 cell_segmentation=Scoreh(:,:,ii)>0;
%                 cell_segmentation=imfill(cell_segmentation,'holes'); %fill holes
%                 cell_segmentation=retain_largest_obj(cell_segmentation); %remove loose parts
%                 cell_segmentation=imclose(imopen(cell_segmentation,ones(3,3)),ones(3,3)); %remove weird protrusions
%                 new_segmentation=new_segmentation+cell_segmentation.*ii; %multiply the segmentation with cell number, so that its location is marked by its number.
%                 new_score=new_score+Scoreh(:,:,ii).*double(cell_segmentation);
%             end
%         else %if there is no splitting event
%             %##Segmentation multiple times##
%             [new_segmentation,new_score]=segmentation_subroutine_for_seeding(subimage_phase,subimage_segmentation,phase_subtraction_factor,dist_modifier,threshold,nrep,max_area_increase_per_tp);
%             new_segmentation=retain_largest_obj(new_segmentation);
%             new_score=retain_largest_obj((new_score.*double(new_segmentation>0)));
%         end
%         %add cell to the Scores and Scores2.
%         final_cell_score=zeros(x_size,y_size,'uint8'); %create an empty matrix that has the same dimensions as the padded phase image.
%         final_cell_score(y_cn,x_cn)=uint8(new_score); %add the cell score
%         Scores(:,:,index)=Scores(:,:,index)+final_cell_score; %##Individual Cell Scores##
%         final_cell_segmentation=zeros(x_size,y_size,'uint8'); %create an empty matrix that has the same dimensions as the padded phase image.
%         final_cell_segmentation(y_cn,x_cn)=uint8(new_segmentation); %add the cell segmentation
%         Scores2(:,:,index)=Scores2(:,:,index)+final_cell_segmentation; %Cell segmentations
%     end %end of cell loop
%     
%     %##The code below corresponds to the pink box in Figure 1D. Above we got
%     %'Individual Cell Scores' held in the matrices Scores and Scores2. Now we will check for overlaps between them and
%     %either merge the cells or distribute the overlaps to get the final results.##
%     
%     %Move the new cells that arose from dividing putative cells to their own Scores sheet.
%     multcellsc=find(numcells>1); %scores with multiple cells
%     for ii=multcellsc
%         for kk=2:numcells(ii)
%             Scoreh=Scores(:,:,ii); %current score
%             Scoreadd=Scores(:,:,ii).*uint8((Scores2(:,:,ii)==kk));
%             Scoreh((Scores2(:,:,ii)==kk))=0; %take it away from the current score
%             Scores=cat(3,Scores,Scoreadd);
%             Scores(:,:,ii)=Scoreh;
%         end
%     end
%     Scores2=(Scores>0); %update cell locations
%     Intersections=sum(Scores2,3); %sum all the 'layers'. If the pixel value is n, n cells are intersecting.
%     Intersections(Intersections==1)=0; %Set the pixels where cells are not intersecting to 0.
%     %go over interactions one by one
%     CC=bwconncomp(Intersections,4); %label individual intersections
%     cells_to_remove=[];
%     cells_to_run=[];
%     %##Check overlaps between cells##
%     for kk=1:CC.NumObjects %go over intersections
%         %find interacting cells=======================
%         Int=zeros([x_size,y_size],'uint8');
%         Int(CC.PixelIdxList{kk})=1; %mask of the kkth intersection
%         Scores_int=Scores.*Int; %multiplies every layer with the kkth intersection.
%         [~,indices]=find(Scores_int); %get the column numbers of the nonzero values
%         indices=unique(indices);
%         [~,j]=ind2sub(y_size,indices); %convert linear indexing to subscripts
%         intersectingcells=unique(j);
%         %===================================================
%         %Merging if necessary
%         for v=1:length(intersectingcells)
%             if ismember(intersectingcells(v),cells_to_remove)==0 %Proceed if the first cell was not removed before.
%                 for w=v+1:length(intersectingcells)
%                     if ismember(intersectingcells(w),cells_to_remove)==0 %Proceed if the second cell was not removed before.
%                         first_cell_segmentation= Scores(:,:,intersectingcells(v))>0;
%                         second_cell_segmentation=Scores(:,:,intersectingcells(w))>0;
%                         Inth= first_cell_segmentation.*second_cell_segmentation; %intersection of the two cells
%                         Intsize=sum(Inth(:)); %size of the intersection
%                         first_cell_size=sum( first_cell_segmentation(:));
%                         second_cell_size=sum(second_cell_segmentation(:));
%                         if max(Intsize/first_cell_size, Intsize/second_cell_size)>0.25 %##Merge?## %If the intersection size is larger than 25% of the smaller cell size then we merge the two cells.
%                             %merge cells by adding the smaller cell to the larger cell and removing the smaller cell.
%                             if first_cell_size>second_cell_size
%                                 %remove the second cell and add parts of the second cell that are not in the first cell to the first cell
%                                 Scores(:,:,intersectingcells(v))=Scores(:,:,intersectingcells(v))+uint8((Scores(:,:,intersectingcells(w)).*uint8(~Inth)));
%                                 cells_to_remove=[cells_to_remove intersectingcells(w)];
%                                 cells_to_run=[cells_to_run intersectingcells(v)]; %run the merged cell again
%                             else
%                                 %the reverse
%                                 Scores(:,:,intersectingcells(w))=Scores(:,:,intersectingcells(w))+uint8((Scores(:,:,intersectingcells(v)).*uint8(~Inth)));
%                                 cells_to_remove=[cells_to_remove intersectingcells(v)];
%                                 cells_to_run=[cells_to_run intersectingcells(w)];
%                                 
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%         
%         
%     end
%     
%     %Running the merged cells again, i.e. feeding these cells into the segmentation_subroutine_for_seeding again to refine the cell boundaries.
%     if isempty(cells_to_run)==0
%         for jj=unique(cells_to_run)
%             no_cell_of_interest=jj;
%             stats=regionprops(logical(Scores(:,:,no_cell_of_interest)),'Centroid','BoundingBox');
%             bbox=round(stats(1).BoundingBox);
%             %The following four lines make sure that the column/row numbers do not exceed the image size.
%             lower_x_limit=max(1,bbox(1)-cell_margin);
%             upper_x_limit=min(y_size,bbox(1)+bbox(3)+cell_margin);
%             lower_y_limit=max(1,bbox(2)-cell_margin);
%             upper_y_limit=min(x_size,bbox(2)+bbox(4)+cell_margin);
%             x_cn=lower_x_limit:upper_x_limit; %column numbers that make up the subimage
%             y_cn=lower_y_limit:upper_y_limit; %row numbers that make up the subimage
%             subimage_phase=double(IPP(y_cn,x_cn));
%             tmpI2=logical(Scores(:,:,no_cell_of_interest));
%             subimage_segmentation=double(tmpI2(y_cn,x_cn));
%             %##Segmentation multiple times## This step refines the cell boundaries.
%             [new_segmentation,new_score]=segmentation_subroutine_for_seeding(subimage_phase,subimage_segmentation,phase_subtraction_factor,dist_modifier,threshold,5,max_area_increase_per_tp);
%             [new_segmentation]=retain_largest_obj(new_segmentation);
%             [new_score]=retain_largest_obj((new_score.*double(new_segmentation>0)));
%             final_cell_score=zeros(x_size,y_size,'uint8');
%             final_cell_score(y_cn,x_cn)=uint8(new_score);
%             Scores(:,:,no_cell_of_interest)=Scores(:,:,no_cell_of_interest)+ final_cell_score;
%         end
%     end
%     
%     Scores(:,:,cells_to_remove)=[];
%     %##Distribute overlaps using scores##
%     [~,Cell_Numbers] = max(Scores,[],3); %compare the score of every cell at every pixel and return the cell number with the highest score at each pixel.
%     %Note that the above line assigns 1 to the pixels where the first cell has the maximum score, but also to the background. To correct that we use the next two lines, which assign zero to background pixels.
%     BG=(Scores(:,:,1)==0).*(~(Cell_Numbers>1));
%     Cell_Numbers(logical(BG))=0; %label matrix, where the location of each cell is marked by its number.
%     %fill the holes
%     number_of_cells=max(Cell_Numbers(:));
%     for i=1:number_of_cells
%         final_cell_segmentation=imfill(Cell_Numbers==i,'holes');
%         final_cell_segmentation=retain_largest_obj(final_cell_segmentation); %in case there are some unconnected pixels.
%         final_cell_segmentation=imclose(imopen(final_cell_segmentation,ones(3,3)),ones(3,3));
%         Cell_Numbers(Cell_Numbers==i)=0;
%         if sum(final_cell_segmentation(:))>min_cell_size %don't include small cells
%             if sum(sum(Cell_Numbers.*double(final_cell_segmentation)))==0
%                 Cell_Numbers=Cell_Numbers+double(i.*final_cell_segmentation);
%             else
%                 Cell_Numbers(logical(Cell_Numbers.*double(final_cell_segmentation)))=0;
%                 Cell_Numbers=Cell_Numbers+double(i.*final_cell_segmentation);
%             end
%         end
%     end
%     %Cell_Numbers is the label matrix, where the location of each cell is marked by its number.
%     %Remove the padding
%     Cell_Numbers(1:size_padding,:)=[]; Cell_Numbers(:,1:size_padding)=[]; Cell_Numbers(end-size_padding+1:end,:)=[]; Cell_Numbers(:,end-size_padding+1:end)=[];
%     %Renumber the cells to account for removed and divided cells.
%     unique_numbers=unique(Cell_Numbers(:));
%     unique_numbers=unique_numbers(2:end)'; %remove 0
%     counter=1;
%     Seed_final=zeros(size(Cell_Numbers));
%     if isempty (unique_numbers)==0
%         for jj=unique_numbers
%             Seed_final=Seed_final+counter.*(Cell_Numbers==jj);
%             counter=counter+1;
%         end
%         
%         
%         %Save the seed.
%         save(save_name,'Seed_final');
%     end
%     %If the code will be run sequentially, comment out the following line.
%     
%     %Visualize the result
%     Cont=zeros(size(Cell_Numbers));%create a matrix to hold the cell contours
%     for i=1:number_of_cells
%         Cont=Cont+bwmorph(Cell_Numbers==i,'remove');
%     end
%     IP_segmentation=uint8(IB); %Create a matrix that will have the segmentation result on the phase image.
%     IP_segmentation(logical(Cont))=255; %Mark the cell contours with the brightest pixel.
%     %figure; imagesc(IP_segmentation)
%     baseFileName = sprintf('Imagefinal #%d.png', timepoint);
%     fullFileName = fullfile('data\output\segmentedcells', baseFileName);
%     imwrite (IP_segmentation, fullFileName)
end