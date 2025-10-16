function batchPlotCorneaPrecipitates(parentdir)
% plots precipitates from all processed folders with xml exports of cornea-scans

list = dir([parentdir filesep '**/*.*']);
subdir = list(cell2mat({list.isdir}));
useful = ~contains({subdir.name},'.');
subdir = subdir(useful);

n = 0;

for i = 1:numel(subdir)
    folder = fullfile(subdir(i).folder, subdir(i).name);
    xml_files = dir(fullfile(folder,'*.xml'));
    precipitatesVolumefile = fullfile(folder,'precipitatesVolume.mat');
    corneafile = fullfile(folder,'cornea.mat');    
    alreadyProcessed = exist(precipitatesVolumefile,'file') & exist(corneafile,'file');

    % if has already been preocessed (i.e. results files in the folder)
    % plot the results
    if alreadyProcessed

        n = n+1;

        % generate an ID from the path, assuming it does not go beyond sub
        % sub dir (2 levels)
        currentParent = subdir(i).folder;
        if length(subdir(i).folder) > length(parentdir)
            ID = [currentParent(length(parentdir)+2:end) '_' subdir(i).name] ;
        else
            ID = [subdir(i).name] ;
        end
        ID = strrep(ID,'\','_');
        ID = strrep(ID,'/','_');
        ID = strrep(ID,', ','_');
        ID = strrep(ID,'(','_');
        ID = strrep(ID,')','_');
        ID = strrep(ID,' ','_');
        ID = strrep(ID,'__','_');
        
        disp(['Showing results for ID: ' ID ]);
        
        load(precipitatesVolumefile);
        load(corneafile);

        xmlfiles = dir(fullfile(folder, '*.xml'));
        S = readstruct( fullfile(folder, xmlfiles(1).name) , 'Filetype', 'xml');
        imagesTags = S.BODY.Patient.Study.Series.Image;
        N = numel(imagesTags)-1;
        w = imagesTags(2).OphthalmicAcquisitionContext.Width;
        h = imagesTags(2).OphthalmicAcquisitionContext.Height;
        x_span = abs(imagesTags(2).OphthalmicAcquisitionContext.End.Coord.X-imagesTags(2).OphthalmicAcquisitionContext.Start.Coord.X);
        y_span = abs(imagesTags(2).OphthalmicAcquisitionContext.Start.Coord.Y-imagesTags(N+1).OphthalmicAcquisitionContext.Start.Coord.Y);
        
        % b-scan (2D) info
        scaleX_2D = imagesTags(2).OphthalmicAcquisitionContext.ScaleX;
        scaleY_2D = imagesTags(2).OphthalmicAcquisitionContext.ScaleY; % this would be scaleZ in 3D (height of the scan)
        pxpmm_x2D = w/x_span; % assuming that the scan is always 11.1 mm wide - this is consitent with the 0.0108 mm/px of scaleX and scaleY
        pxpmm_y2D = 1/scaleY_2D;
        
        % volume (3D) info
        scaleX_3D = scaleX_2D;
        scaleY_3D = y_span/(N-1);
        scaleZ_3D = scaleY_2D;
        pxpmm_x3D = pxpmm_x2D;
        pxpmm_y3D = 1/scaleY_3D;
        pxpmm_z3D = pxpmm_y2D;
        
        scanStart1 = [imagesTags(2).OphthalmicAcquisitionContext.Start.Coord.X      imagesTags(2).OphthalmicAcquisitionContext.Start.Coord.Y];
        scanEnd1 = [imagesTags(2).OphthalmicAcquisitionContext.End.Coord.X          imagesTags(2).OphthalmicAcquisitionContext.End.Coord.Y];
        scanStartN = [imagesTags(N+1).OphthalmicAcquisitionContext.Start.Coord.X    imagesTags(N+1).OphthalmicAcquisitionContext.Start.Coord.Y];
        scanEndN = [imagesTags(N+1).OphthalmicAcquisitionContext.End.Coord.X       imagesTags(N+1).OphthalmicAcquisitionContext.End.Coord.Y];
        
        fixedPoints = ([ scanStartN*pxpmm_x2D ; scanEndN*pxpmm_x2D;  scanStart1*pxpmm_x2D ; scanEnd1*pxpmm_x2D  ]);
        movingPoints = fliplr([N  1; N w; 1 1; 1 w]);
        tform = fitgeotrans(movingPoints,fixedPoints,'affine');
        
        movingPointsIntrp = fliplr([round(y_span/11.1*1024)  1; round(y_span/11.1*1024) w; 1 1; 1 w]);
        tformIntrp = fitgeotrans(movingPointsIntrp,fixedPoints,'affine');
        
        distanceFromTop = .5; % relative to scan height
        
        % localizer image
        [~,name,ext] = fileparts(imagesTags(1).ImageData.ExamURL);
        imfile = fullfile(folder, [char(name) char(ext)]);
        localizer = rgb2gray( imread(imfile) );
                
        %% fitted boundary of the analysed corneaInt area (1/2 from top of the scan) 
        slice = abs(corneaInt-round(h*distanceFromTop));
        [~, min1] = min(slice,[],2);
        slice(1:N, min1) = h;
        [~, min2] = min(slice,[],2);
        lh = size(localizer,1);
        loc_y_bot = scanStart1(2);
        loc_y_top = scanStartN(2);
        loc_y_botpx = loc_y_bot/11.1*lh;
        loc_y_toppx = loc_y_top/11.1*lh;
        range1toNinPx = fliplr(loc_y_toppx : (loc_y_botpx-loc_y_toppx)/(N-1) : loc_y_botpx);
        slcX= [min1; min2];
        slcY = [range1toNinPx range1toNinPx]';
        % meanSlc = mean(slcX);
        % slcL = [slcX(slcX<meanSlc) slcY(slcX<meanSlc)];
        % slcR = [slcX(slcX>meanSlc) slcY(slcX>meanSlc)];
        % unfinished ... 
        
        %% surf plot 
        % prec = precipitates;
        % prec(isnan(prec)) = 0;
        % figure; s = surf(242-corneaInt-prec);
        % s.EdgeColor = 'none';
        % colormap jet;
        % s.FaceColor = 'interp';
        % s.FaceAlpha = .5;
        % hold on;
        % scatter3(slcX, [1:81 1:81]', repmat(121, [numel(slcX) 1]),'filled', 'MarkerFaceColor', [0 1 0]);
        % 
        % xlim([0 1024]);
        % xticks((0:11).*pxpmm_x3D);
        % xticklabels( {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'});
        % xlabel('mm');
        % 
        % yticks((0:5).*pxpmm_y3D);
        % yticklabels( {'0', '1', '2', '3', '4', '5'});
        % ylabel('mm');
        % 
        % zlim([0 242]);
        % zticks((0:2).*pxpmm_z3D);
        % zticklabels( {'0', '1', '2'});
        % zlabel('mm');
        
        %% show precipitates on localizer with roi boundary with the un-interpolated precipiates map
        
        % figure; imshow(gray2rgb(localizer));
        % hold on
        % canvas = imref2d(size(localizer));
        % warpedPrec = imwarp(precipitates, tform, 'OutputView', canvas );
        % OverlayImage = imshow( warpedPrec );
        % caxis auto
        % colormap( OverlayImage.Parent, jet );
        % colorbar( OverlayImage.Parent);
        % cb = colorbar( OverlayImage.Parent);
        % cb.Ticks = [0 0.05 0.1 0.15]./scaleZ_3d;
        % cb.TickLabels = ["0" "0.05" "0.1" "0.15"];
        % cb.Label.String = "Height (mm)";
        % set(OverlayImage, 'AlphaData', double( warpedPrec./10 ));
        % plot(slcX,slcY,'g.');
        % volumeBound = [scanStart1; scanEnd1; scanEndN; scanStartN; scanStart1]*pxpmm_x2D;
        % volumeBound(volumeBound==0) = 1;
        % plot(volumeBound(:,1), volumeBound(:,2), 'g-', 'LineWidth', 2);
        
        %% show precipitates on localizer with roi boundary with the interpolated precipiates map
        distanceFromTop = .5; % relative to scan height
        prec = precipitates;
        analysedRegion = (corneaInt<h*distanceFromTop);
        prec(isnan(prec)&analysedRegion) = 0;
        precClean = prec;
        % cleaning the precipitates' map from smaller detections 
        % precClean(precClean<1) = 0;
        % cleaning the precipitates' map from reflection-artifacts 
        % precClean(precClean>30) = 0;
        
        % interpolating the map
        precIntrp = imresize(precClean, [round(size(precClean,1)*scaleY_3D/scaleY_2D) size(precClean,2)], 'bicubic' ).*scaleZ_3D;
        precIntrpPlot = precIntrp;
        precIntrp(precIntrp==0) = nan;
        
        % showing the map on the localiser
        fh = figure;
        imshow(gray2rgb(localizer));
        hold on
        canvas = imref2d(size(localizer));
        warpedPrec = imwarp(precIntrpPlot, tformIntrp, 'OutputView', canvas );
        warpedPrec(1:round(fixedPoints(1,2)),:) = nan;
        warpedPrec(round(fixedPoints(3,2)):end,:) = nan;
        OverlayImage = imshow( warpedPrec );
        caxis auto
        colormap( OverlayImage.Parent, jet );
        colorbar( OverlayImage.Parent);
        cb = colorbar( OverlayImage.Parent);
        
        % colormap to be used if precipitates' map has not yet been scaled to
        % % physical scale (scaleZ_3D)
        % cb.Ticks = [0 0.05 0.1 0.15]./scaleZ_3D;
        % cb.TickLabels = ["0" "0.05" "0.1" "0.15"];
        % set(OverlayImage, 'AlphaData', double( warpedPrec./10 ));
        
        % colormap to be used if precipitates' map has been scaled to physical
        % scale (scaleZ_3D)
        zero_off = .02; % offset of transparency at zero 
        caxis([0 .14])
        cb.Ticks = [0 : 0.02 : 0.14];
        set(OverlayImage, 'AlphaData', double( (warpedPrec+zero_off).*10 ));
        cb.Label.String = "Elevation (mm)";
        cb.LineWidth = 1;
        cb.FontSize = 12;
        
        % plot(slcX,slcY,'g.');
        volumeBound = [scanStart1; scanEnd1; scanEndN; scanStartN; scanStart1]*pxpmm_x2D;
        volumeBound(volumeBound==0) = 1;
        plot(volumeBound(:,1), volumeBound(:,2), 'g-', 'LineWidth', 1);
        set(gcf, 'units','normalized', 'position',[.5 .1 .5 .7])
        set(gca, 'units','normalized', 'innerposition',[0 .05 1 .9])
        % To save or print a figure with specific dimensions, set the 
        % PaperPosition property of the figure to the desired dimensions. 
        % The PaperPosition property affects the size of saved and printed 
        % figures, but does not affect the size of the figure on the 
        % display. Set the property to a four-element vector of the form 
        % [left bottom width height].
%         fh.InvertHardcopy = 'off'; % Ensure that the colors of the saved figure match the colors on the display by setting the InvertHardcopy property of the figure to 'off'. 
%         saveas(gcf, fullfile(parentdir, [ID '.png']) );
        print(fh,fullfile('C:\Users\Giovanni\Dropbox\Precipitati cheratici corneali\AS-OCT for gio\XML\all_result_images', [ID '.png']),'-dpng','-r300');
    end
end



