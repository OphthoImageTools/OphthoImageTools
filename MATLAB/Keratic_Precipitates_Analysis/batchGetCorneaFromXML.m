function [corneaExt, corneaInt, precipitates, enface, images]  = batchGetCorneaFromXML(folderWithXml)

showScansWhileProcesing = true;
noXmlError = 'No .xml file found in the given folder. Cannot continue, operations terminated.';
multipleXmlError = 'The folder contains more than one .xml file. Cannot continue, operations terminated.';

xmlfiles = dir(fullfile(folderWithXml, '*.xml'));

if numel(xmlfiles) > 1 
    error(multipleXmlError);
elseif numel(xmlfiles) == 0 
    error(noXmlError)
end

S = readstruct( fullfile( folderWithXml, xmlfiles(1).name) , 'Filetype', 'xml');

% read the info about the scan size
imagesTags = S.BODY.Patient.Study.Series.Image;
N = numel(imagesTags)-1;
w = imagesTags(2).OphthalmicAcquisitionContext.Width;
h = imagesTags(2).OphthalmicAcquisitionContext.Height;
scaleX = imagesTags(2).OphthalmicAcquisitionContext.ScaleX;
scaleY = imagesTags(2).OphthalmicAcquisitionContext.ScaleY;
pxpmm = 1024/11.1;

% warning off to keep console clean of fitting warnings (very likely)
warning off

f = waitbar(0,'Please wait...');

cint = [];
prec = [];
corneaInt = nan(N,w);
corneaExt = nan(N,w);
precipitates = nan(N,w);
precipitates0 = nan(N,w);
enfaceDist = 20; % pixels
enface = nan(N,w);

tic 
for i =  1:N+1%round(N*.92):N+1

    waitbar(i/(N+1), f, ['processing : ' num2str(i) '/' num2str(N+1)] );
    

    try
        
        % get the file
        [~,name,ext] = fileparts(imagesTags(i).ImageData.ExamURL);
        imfile = fullfile(folderWithXml, [char(name) char(ext)]);

        if exist(imfile, "file")
            if strcmpi(char(imagesTags(i).ImageType.Type), "LOCALIZER")
                disp([ num2str(i-1) ': ' char(name) char(ext) ' - LOCALIZER']);
    
                % save the infrared image (localizer)
                imIR = rgb2gray( imread(imfile) );
            
            elseif strcmpi(imagesTags(i).ImageType.Type, "OCT")
    
                disp([ num2str(i-1) ': ' char(name) char(ext)]);
                
                % save scan position relative to the localizer
                scanStart{i-1} = [imagesTags(i).OphthalmicAcquisitionContext.Start.Coord.X imagesTags(i).OphthalmicAcquisitionContext.Start.Coord.Y];
                scanEnd{i-1} = [imagesTags(i).OphthalmicAcquisitionContext.End.Coord.X imagesTags(i).OphthalmicAcquisitionContext.End.Coord.Y];
    
                % read and process the scans
                im{i-1} =  rgb2gray( imread(imfile) );
                if size(im{i-1},1) == 242
                    images(:,:,i-1) = im{i-1} ; 
                    [cext, cextf, mask, noiseMean, noiseStd]  = getCorneaExternalBoundary(im{i-1});
                    if ~isempty(cextf) && ~isempty(noiseMean) && ~isempty(noiseStd)
                        [cint, prec, mask] = getCorneaInternalBoundary(im{i-1}, cextf,  mask, noiseMean, noiseStd);
                    end
        
                    % take the ext boundary if found
                    if ~isempty(cext)
                        corneaExt(i-1, cext(:,1)) = cext(:,2);
        
                    % or snake using as guess the previous boundary, if found 
                    elseif ( i>2 && any(~isnan(corneaExt(i-2,:))) ) 
                        disp('Trying to use the boundary from the previous scan.')
        
                        % snake 1D from the fitted curve
                        Options.Wline=.01;                % 0
                        Options.Wedge=4;               % 10
                        Options.Sigma1=round(2);  % 2
                        Options.Sigma2=round(2);  % 3
                        Options.Alpha=1;               % .6
                        Options.Beta=10;              % 10^4
                        Options.Gamma=1;                % 1
                        Options.Kappa=10;               % .4
                        Options.Iterations=100;          % 20
                        Options.B2Wtransition = 1;
                        Options.Verbose = false;
                        corneaExtguess = corneaExt(i-2, :)-10 ;
                        corneaExtguess(isnan(corneaExtguess)) = size(images,1);
                        corneaExtguess = [corneaExtguess' (1:length(corneaExtguess))'];
                        corneaExtSnk = Snake1D(im{i-1},corneaExtguess,Options);
        
                        mask = true (size(images,1), size(images,2));
                        notMask = magicwandmono( imadjust(images(:,:,i-1),[],[],.2), [1 1 size(images,1) size(images,1)], [1 size(images,2) size(images,2) 1], .02);
                        mask(notMask) = false;
                        closedMask = imfill(imclose(mask,strel('disk',3)),'holes');
                        mask = imopen(closedMask,strel('disk',3));
                        x = corneaExtSnk(1,:);
                        y = corneaExtSnk(2,:);
                        inMask = pointsInsideMask(mask,x,y);
                        cext = [x(inMask) y(inMask)];
            
                        corneaExt(i-1, cext(:,1)) = cext(:,2);
                    end
        
                    % take the int boundary if found
                    if ~isempty(cint)
                        corneaInt(i-1, cint(:,1)) = cint(:,2);
                        enface(i-1,:) = sum( im{i-1} .* uint8(maskBelow(im{i-1},corneaInt(i-1,:) ) & maskAbove(im{i-1},corneaInt(i-1,:)+20 )), 1);  
                    % or snake using as guess the previous boundary, if found 
                    elseif ( i>2 && any(~isnan(corneaInt(i-2,:))) ) 
                        disp('Trying to use the boundary from the previous scan.')
                        xdata = find(corneaInt(i-2,:)>0)';
                        ydata = corneaInt(i-2,xdata)';
                        f_prev = fit(xdata,ydata,'poly3');
                        corneaIntf = feval(f_prev,1:size(im{i-1},2)); 
                        [cint, prec, mask] = getCorneaInternalBoundary(im{i-1}, cextf-40, mask, noiseMean, noiseStd, corneaIntf);
                        if ~isempty(cint)
                            disp("Found with the boundary of the preivouse scan.")
                            corneaInt(i-1, cint(:,1)) = cint(:,2);
                            enface(i-1,:) = sum( im{i-1} .* uint8(maskBelow(im{i-1},corneaInt(i-1,:) ) & maskAbove(im{i-1},corneaInt(i-1,:)+20 )), 1); 
                        else
                            disp("Still unable to find it.")
                        end
                    end
        
                    if ~isempty(prec)
                        precipitates0(i-1,:) = sum(prec,1);
                        p0 = precipitates0(i-1,:);
                        precipitates(i-1,p0>0) = p0(p0>0);
                    end
    
                else
                    disp('Error on the height of the scan, skipping');
                end
    
                % only show results if showScansWhileProcesing == true
                if showScansWhileProcesing
                    figure(1)
                    imshow(im{i-1});
                    hold on
                    if ~isempty(cext)
                        plot(corneaExt(i-1, :),'r-');
                    end
                    if ~isempty(cint)
                        plot(corneaInt(i-1, :),'b-');
                    end
                    if ~isempty(prec)
                        precOnBndry = corneaInt(i-1, :)+precipitates0(i-1,:);
                        precOnBndry(p0<=0) =nan;
                        plot(precOnBndry,'g-');
                    end
                    plot(1:size(images,2), repmat(size(images,1)/2,1,size(images,2)), '-m');
    %                 ymean = mean(im{i-1},1); 
    %                 aboveAVG = ymean > mean(ymean)+2*std(ymean); 
    %                 plot(ymean);
    %                 plot(find(aboveAVG), ymean(aboveAVG));
                    hold off
                end
                
            end
        else
            disp([imfile  ' not found'] );
        end

    catch
        disp(['Error processing ' num2str(i-1) ': ' char(name) char(ext) ])
    end
    
end
disp(toc)


% FOLLOWING OPERATION IS NOT RELIABLE. BETTER TO EXCLUDE SCANS WHERE IT WAS
% NOT FOUND, AS IT IS OFTEN ASSOCIATED WITH ARTIFACTS/BAD SCANS AND
% MEASUREMENTS ARE NOT RELIABLE. IT IS BETTER TO IGNORE THESE SCANS AND
% MITIGATE THEIR EFFECT IN THE MEASUREMENTS DIVIDING BY THE TOTAL
% SURFACE AREA OF INTERNAL BOUNDARY - WHERE THESE SCANS ARE EXCLUDED
% get internal boundary with linear interp where not found after the first
% scan where it was found
corneaIntNotFound = ~any(~isnan(corneaInt),2);
firstScanWithCorneaInt = find(~corneaIntNotFound,1);
corneaIntNotFound(1:firstScanWithCorneaInt-1)=false;

% if any(corneaIntNotFound)
%     corneaIntFilled =  fillmissing(corneaInt,'linear');
%     corneaInt(corneaIntNotFound,:) = corneaIntFilled(corneaIntNotFound,:);
%     
%     % find the precipitates for those scans without found cInt that has been
%     % filled with "fillmissing"
%     for k = find(corneaIntNotFound)'
%         disp(['calculating precipitates for scan ' num2str(k)])
%         scan = im{k};
%     
%         mask = true (size(scan,1), size(scan,2));
%         notMask = magicwandmono( imadjust(scan,[],[],.2), [1 1 size(scan,1) size(scan,1)], [1 size(scan,2) size(scan,2) 1], .02);
%         mask(notMask) = false;
%         closedMask = imfill(imclose(mask,strel('disk',3)),'holes');
%         mask = imopen(closedMask,strel('disk',3));
%         erodedMask = imerode(mask,strel('disk',10));
%     
%         cext_k = corneaExt(k,:);
%         cext_k(isnan(cext_k)) = size(scan,1);
%     
%         % calculate the noise intensity
%         air = scan(maskAbove(scan,cext_k) & mask);
%         noiseMean = mean(air);
%         noiseStd = round(std(double(air)));
%     
%         % cornea-enhanced image
%         imeq = imadjust(scan,[0 (noiseMean+2*noiseStd)/255 ],[0 1]);
%         corneaEnh_im = gaussianbpf(imeq,2,50);
%         corneaEnh_im = corneaEnh_im == 255;
%         corneaEnh_im = imfill(corneaEnh_im,8,'holes');
%     
%         % approach based on canny edge detection
%         BW =  edge(corneaEnh_im,'canny',.8, 3) ;
%         BWmasked = BW ;
%     
%         % take only edges under the external boundary and above the bottom edge of
%         % the scan 
%         BWmasked(~maskBelow(mask, cext_k) | ~erodedMask) = 0;
%         BWmasked = bwskel(BWmasked);
%         
%         % find the biggest component of the edges to find the top of the internal
%         % boundary of the cornea, to calculate the search area. The search area is
%         % calculated from the external boundary, based on the shape and distance of
%         % the external boundary (with increased curve, as calculated below).
%         CC = bwconncomp(BWmasked);
%         numPixels = cellfun(@numel,CC.PixelIdxList);
%         [~,idx] = max(numPixels);
%     
%         % % take only what is below the max component
%         [ydata, xdata] = ind2sub(size(BW), CC.PixelIdxList{idx});
%         [topInternal, locMinInt] = min(ydata);
%         [topExternal, locMinExt] = min(cext_k);
%         distIntExt =  topInternal - topExternal;
%         curveIncrease = max( 2, (1+abs(distIntExt-55)/5) )  * 100 * ( 1-sin(acos(([1:1024] - round((locMinExt+xdata(locMinInt))/2 ) )   / 1024)));
%         maskInt = ( maskBelow(mask, cext_k+distIntExt-5+curveIncrease) & maskAbove(mask, 40+distIntExt+cext_k+curveIncrease) & erodedMask);
%         BWmasked = maskInt & BW;
%     
%         % from the search area take only the first (column-wise) edge from the bottom of the scan
%         [~, loc] = max( flipud(BWmasked), [], 1 );
%         cannyInt = 1+size(im,1)-loc;
%         cannyInt(cannyInt >= size(im,1)) = 0;
% 
%         % find precipitates as the areas above the canny edges and below the
%         % internal boundary of the cornea ( and above an arbitrary offset -
%         % should this be half of the scan?)
%         cint_k = corneaInt(k,:);
%         cint_k(isnan(cint_k)) = size(scan,1);
%         prec_k = maskBelow(mask, cint_k) & maskAbove(mask,cannyInt) &  [true([size(scan,1)/2 size(scan,2)]); false([size(scan,1)/2 size(scan,2)])] ;
%         precipitates0(k,:) = sum(prec_k,1);
%         p0 = precipitates0(k,:);
%         precipitates(k,p0>0) = p0(p0>0);
%     end
% end

close(f);
warning on;

% [X,Y] = meshgrid( 1:w, 1:N);
% figure;
% % surf cornea int
% s = surf(X,Y,h-corneaInt,'FaceAlpha',0.5);
% s.EdgeColor = 'none';
% % surf cornea int with precipitates
% sp = surf(X,Y,h-corneaInt-precipitates0,'FaceAlpha',0.5);
% sp.EdgeColor = 'none';
% 
% 
% 
% enface_im = (uint8(255*enface/max(enface(:))));
% 
% fixedPoints = ([ scanStart{N}*pxpmm ; scanEnd{N}*pxpmm;  scanStart{1}*pxpmm ; scanEnd{1}*pxpmm  ]);
% movingPoints = fliplr([N  1; N w; 1 1; 1 w]);
% tform = fitgeotrans(movingPoints,fixedPoints,'affine');
% Roriginal = imref2d(size(imIR));
% enface_im = imwarp(enface_im, tform, 'OutputView',Roriginal);
% figure; imshow(imIR);
% hold on; 
% r = imshow(enface_im);
% set(r, 'AlphaData', enface_im); 
% 
% 
% precipitates_im = imwarp(precipitates0, tform, 'OutputView',Roriginal);
% figure; imshow(imIR);
% hold on; 
% r = imshow(precipitates_im);
% set(r, 'AlphaData', precipitates_im); 
% 
% 
% figure; p = surf(X,Y,precipitates0,'FaceAlpha',0.5);
% p.EdgeColor = 'none';
% p.FaceColor = 'interp';
% zlim ([0 50])


% zLocs = 1:(w/2-1)/(N-1):w/2;
% 
% x = [];
% y = [];
% z = [];
% for i = 1 : N
%     prec = imdilate(imerode(precipitates{i},strel('disk',1)),strel('disk',1));
%     [row,col] = find(prec);
%     y = [y; row];
%     x = [x; col];
%     z = [z; repmat(zLocs(i), size(row))];
% end
% 
% corneaExt3D = [];
% for i = 1 : N
%     corneaExt3D = [corneaExt3D; [corneaExt{i} repmat(zLocs(i), [size(corneaExt{i},1) 1] )]];
% end
% 
% surffit = fit([corneaExt3D(:,1),corneaExt3D(:,3)],corneaExt3D(:,2),'poly22','normalize','on');
% [X,Y] = meshgrid(1:(w-1)/(N-1)*2:w, 1:(w/2-1)/(N-1):(w/2));
% Z = feval(surffit,X,Y);
% figure;
% s = surf(X,Y,Z,'FaceAlpha',0.3);
% s.EdgeColor = 'none';
% axis equal
% colormap winter
% hold on
% sc = scatter3(x,512-z,y,'filled');
% sc.SizeData = 5;