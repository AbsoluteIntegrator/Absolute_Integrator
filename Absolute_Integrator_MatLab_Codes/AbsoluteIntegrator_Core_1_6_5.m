% Close all waitbars.
F = findall(0,'type','figure','tag','TMWWaitbar'); % Find all waitbars...
delete(F);
clear F

if exist('now','var') == 1  % This line prevens circumvention of date check by the user specifying a variable called 'now'.
    clear now
end

run Integrator_license_check_1_0.m

% Check if the date IS beyond the expiry AND the user IS offline:
if floor(now) > 736542 && strcmp(licenseStatus,'Online - valid') == 0 % Date should match here and in warning string below.
    % Display expiration text:
    msgbox (horzcat('This offline copy of ',version_string,' has expired. Please connect to the internet and retry. For a new license please contact support@lewysjones.com .'), version_string , 'warn')
else
    % If either of the conditions above are false execute the main code:
    
    % If the user is offline, display warning message about expiration:
    if strcmp(licenseStatus,'Online - valid') == 0
        expiration_string = datestr(736542) ; % Date should match here and in logic check above.
        expiration_string = horzcat('This offline copy of ',version_string,' will expire on ',expiration_string,'. To upgrade your license please contact support@lewysjones.com .') ;
        h = msgbox ( expiration_string , version_string , 'help' );
        edithandle = findobj(h,'Style','pushbutton') ;
        set(edithandle,'Visible','off') ;
        pause(3)
        close (h)
        pause(0.05)
        clear expiration_string
        clear edithandle
        clear h
    end
    
    %*********************** MAIN CODE BEGINS HERE ************************
    
    % Variable pre-handling:
    format compact
    currentNormalisation   = faradayCurrentImaging / faradayCurrentDetector ; % Increse in current used for imaging.
    
    if manual_best_size == 0
        clear best_size
        clear manual_best_size
    end
    
    if loadSimulation == 1
        currentNormalisation  = 1 ; % Verification of current setting for simulation.
        low_lim = 0 ;
    end
    
    if use_SpectrumIntegrator == 1 && upsamplingFactor ~= 1
        warning('You have selected to use SpectrumIntegrator, image upsampling cannot be used. Upsampling disabled, you may need to respecify peak-finding settings.')
    end
    
    % END Variable pre-handling.
    
    % Tell user how many files to open and in what sequence:
    if loadSimulation == 1
        number_files = 1 ;
        message = ['Please open ', num2str(number_files) , ' file when prompted:' , char(10) , '--> File should be the simulated image.'] ;
    elseif analyse_CCD_flux == 0
        number_files = 2 ;
        message = ['Please open ', num2str(number_files) , ' files when prompted:' , char(10) , '--> First file should be the detector response map,' , char(10) , '   --> Second file should be the experimental image.'] ;
    else
        number_files = 3 ;
        message = ['Please open ', num2str(number_files) , ' files when prompted:' , char(10) , '-->  First file should be the detector response map,' , char(10) , '   --> Second file should be the experimental CCD Flux image,' , char(10) , '      --> Third file should be the experimental image.'] ;
    end
    
    clear number_files
    uiwait(msgbox(message,version_string));
    
    
    % End variable pre-handling, main code begins.
    
    fullscreen = get(0,'ScreenSize'); % This line determines the size of the users screen:
    
    if loadSimulation == 0
        % User loads detector image:
        h = msgbox('Loading detector image and finding center. May take > 1min for large files.',version_string,'help') ;
        run 'C:\Users\Lewys\Box Sync\DPhil\MatLab Software\File Opener\File_Opener_1_5_5'
        close(h)
        detectorImage = input_file ;
        
        %         warning('Detector cropped here')
        %         clear detectorImage
        %         detectorImage = input_file(:,1:256) * 0.1 ;
        %     min_detector = min(min(detectorImage))
        %     mean_detector = mean(mean(detectorImage))
        
        clear input_file
        
        %         warning('Detector scan resizing disabled!!')
        if size(detectorImage) >= 2048
            detectorImage = imresize(detectorImage, 0.25 );
        elseif size(detectorImage) >= 1024
            detectorImage = imresize(detectorImage, 0.5 );
        elseif size(detectorImage) < 384
            detectorImage = imresize(detectorImage, 2 );
        end
        
        if weightDetector == 1
            
            % First a fractional detector image is created that is thresholded at
            % the 5% level:
            detectorThresholded = (detectorImage - min(min(detectorImage))) / max(max(detectorImage - min(min(detectorImage))));
            detectorThresholded = medfilt2(detectorThresholded) ;
            detectorThresholded(detectorThresholded >= 0.3) = 1 ;
            detectorThresholded(detectorThresholded < 0.3) = 0 ;
            activeLayer = detectorThresholded ;             % Saves the deduced active layer for weighting later.
            detectorThresholded = 1 - detectorThresholded ; % Changed to be one in vaccum to find centre point.
            detectorThresholded = bwmorph(detectorThresholded,'close') ;
            detectorThresholded = bwmorph(detectorThresholded,'clean') ;
            
            SE = strel('disk', round(0.01*size(detectorImage,1)), 0) ;
            black_area = ( 1 - imdilate(activeLayer,SE)) .* detectorImage ; % Adds a small dilation to avoid mask edge clashes.
            clear SE
            
            black_pix = nonzeros(black_area) ;
            clear black_area
            vacuumIntesity = median(black_pix) ;
            clear black_pix
            
            % Image erosion to create points at likely detector centres:
            detectorErroded = bwulterode(detectorThresholded) ;
            clear detectorThresholded
            
            % Fiding of likely detector centres:
            [row,col] = find(detectorErroded) ;
            
            % Ranking of distance from image centre:
            imageCentre = size(detectorErroded) / 2 ;
            clear detectorErroded
            distance = sqrt( (row-imageCentre(1)).^2 + (col-imageCentre(2)).^2 ) ;
            clear imageCentre
            [~,I] = min(distance) ;
            centrePoint = [row(I) col(I)] ;
            clear distance
            clear row
            clear col
            
            
            
            % Next create polar coordinates centered on the detector centre:
            [r, c] = size(detectorImage) ;
            [X, Y] = meshgrid(1:c,1:r) ;
            clear r
            clear c
            X = X - centrePoint(2) ;
            Y = Y - centrePoint(1) ;
            [theta , rho] = cart2pol(X, Y); % Here 'rho' is measuered in pixels but should be converted to milli-radians later.
            theta = theta * 180 / pi ;
            clear X
            clear Y
            
            
            h = waitbar(0,'Analysing Radial Sensitivity. Please wait...','Name',version_string);
            
            rhoMax = round(max(rho(:))) ; % The rho (in pixels) to stop at.
            meanSensitivity = zeros(1,rhoMax) ;
            for angularRange = 0:(rhoMax-1) % Steps through in pixels radius range from zero-one, one-two etc..
                angleMask = ones(size(detectorImage)) ; % Creates mask of ones in chosen range.
                angleMask(rho <= angularRange) = 0 ;
                angleMask(rho > angularRange+1) = 0 ;
                visibleDetector = angleMask .* (detectorImage-vacuumIntesity) ; % Applies mask.
                visibleDetector(visibleDetector == 0 ) = [] ; % Strips zeros.
                meanSensitivity(angularRange+1) = mean(mean(visibleDetector)) ;
                waitbar(angularRange / rhoMax ) ;
            end
            clear angleMask
            
            delete(h)
            
            sensitivityGradient = gradient(meanSensitivity) ;
            [~,I] = max(sensitivityGradient) ; % Finds pixel index of inner angle onset.
            [~,O] = min(sensitivityGradient) ; % Finds pixel index of outer angle.
            clear sensitivityGradient
            mradPerPixel = innerAngle / I ;
            
            if analyseAzimuthalSensitivity == 1
                
                h = waitbar(0,'Analysing Azimuthal Sensitivity. Please wait...','Name',version_string);
                azimuthalSensitivity = zeros(1,361) ;
                for angularRange = -180:180 % Steps through in angle range...
                    angleMask = ones(size(detectorImage)) ; % Creates mask of ones in chosen range.
                    angleMask(theta <= angularRange-0.5) = 0 ;
                    angleMask(theta > angularRange+0.5) = 0 ;
                    angleMask(rho < 1.1*I) = 0 ;
                    angleMask(rho > 0.75*O) = 0 ;
                    visibleDetector = angleMask .* (detectorImage-vacuumIntesity) ; % Applies mask.
                    azimuthalSensitivity(angularRange+181) = mean(mean(visibleDetector)) ;
                    waitbar( (angularRange+181) / 360 ) ;
                end
                
                azimuthalSensitivity = medfilt1(azimuthalSensitivity) ;
                azimuthalSensitivity = azimuthalSensitivity / mean(azimuthalSensitivity(:)) ;
                
                roundnessScore = 100 * std(azimuthalSensitivity)
                
                figure('Position',[108 148 fullscreen(3)-250 fullscreen(4)-396],'Name',horzcat('Detector Azimuthal Sensitivity - ',version_string),'NumberTitle','off' ) ;
                bar(0:360,azimuthalSensitivity,'FaceColor',[.714 0.475 1])
                xlim([0 360])
                set(gca,'FontSize',28,'LineWidth',2,'XColor','k','Color',[0.914 0.914 0.914])
                set(gca,'layer','top')
                set(gca, 'GridLineStyle', '-');
                hline = refline([0 1]);
                set(hline,'Color','r','LineWidth',3)
                xlabel('Azimuthal Angle (degrees)','FontSize',34,'FontWeight','bold')
                ylabel('Relative Sensitivity','FontSize',34,'FontWeight','bold')
                
                clear angleMask
                clear visibleDetector
                clear rhoMax
                clear theta
                delete(h)
                
            end
            
            
            
            detectorAngleAxis = linspace( 0 , (mradPerPixel*(size(meanSensitivity,2))) , size(meanSensitivity,2) ) ;
            
            % Ooption to analyse CCD flux image to determine appropriate TDS exponent:
            if analyse_CCD_flux == 1
                % The flux-image will be analysed and the variable 'TDSexponent' will be overwritten.
                run Analyse_CCD_Flux_0_4
            end
            
            figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Detector Centre Location - ',version_string),'NumberTitle','off' ) ;
            subplot(1,4,1)
            imagesc(detectorImage-vacuumIntesity)
            set(gca,'XTick',[],'YTick',[])
            axis image
            title('Detector Sensitivity (counts)','FontSize',11,'FontWeight','bold')
            colorbar('horiz');
            hold on
            plot(centrePoint(2),centrePoint(1),'p','MarkerEdgeColor','r','MarkerFaceColor','r')
            hold off
            
            subplot(1,4,2)
            imagesc(activeLayer)
            set(gca,'XTick',[],'YTick',[])
            axis image
            title('Detector Active Area','FontSize',11,'FontWeight','bold')
            colorbar('horiz');
            hold on
            plot(centrePoint(2),centrePoint(1),'p','MarkerEdgeColor','r','MarkerFaceColor','r')
            hold off
            
            fluenceCone = (rho*mradPerPixel).^-TDSexponent ;  % Generates the TDS fluence distribution.
            fluenceCone(fluenceCone==Inf) = 0        ;  % Strip out the central pixel.
            fluenceCone = activeLayer .* fluenceCone ;  % Only consideres the fluence falling on active areas of the detector.
            
            
            if analyse_CCD_flux == 1
                fluenceCone( rho*mradPerPixel > flux_cuttof_angle ) = 0        ;  % Strip out the central pixel.
            end
            fluenceTotal = sum(fluenceCone(:))       ;  % Calcualte the total fluence for on axis detector.
            fluenceCone = fluenceCone / fluenceTotal ;  % Normalises the total remaining fluence over the active region to unity.
            
            if sum(sum(offsetDetector.^2)) > 0
                % Recreate polar coordinates that are displaced:
                [r, c] = size(detectorImage) ;
                [X, Y] = meshgrid(1:c,1:r) ;
                clear r
                clear c
                X = X - centrePoint(2) - offsetDetector(2)*I/100 ;
                Y = Y - centrePoint(1) - offsetDetector(1)*I/100 ;
                [~ , rho] = cart2pol(X, Y); % Here 'rho' is measuered in pixels but should be converted to milli-radians later.
                clear X
                clear Y
                
                % Recalcualte the fluence cone:
                fluenceCone = (rho*mradPerPixel).^-TDSexponent ;  % Generates the TDS fluence distribution.
                fluenceCone(fluenceCone==Inf) = 0        ;  % Strip out the central pixel.
                fluenceCone = activeLayer .* fluenceCone ;  % Only consideres the fluence falling on active areas of the detector.
                fluenceCone = fluenceCone / fluenceTotal ;  % Normalises the total remaining fluence usingf the previous esitmate.
            end
            
            clear rho
            clear activeLayer
            
            subplot(1,4,3)
            imagesc(fluenceCone)
            set(gca,'XTick',[],'YTick',[])
            axis image
            title('Electron Flux Distribution','FontSize',11,'FontWeight','bold')
            colorbar('horiz');
            hold on
            plot(centrePoint(2),centrePoint(1),'p','MarkerEdgeColor','r','MarkerFaceColor','r')
            hold off
            
            fluenceWeightedDetector =  fluenceCone .* (detectorImage-vacuumIntesity) ;   % Applies the fluence weighting over the detector.
            clear fluenceCone
            
            subplot(1,4,4)
            imagesc(fluenceWeightedDetector)
            set(gca,'XTick',[],'YTick',[])
            axis image
            title('Flux Weighted Sensitivity','FontSize',11,'FontWeight','bold')
            colorbar('horiz');
            hold on
            plot(centrePoint(2),centrePoint(1),'p','MarkerEdgeColor','r','MarkerFaceColor','r')
            hold off
            
            beamIntensity = sum(fluenceWeightedDetector(:)) ;   % Integrates the sensitivity to find the fluence weighted mean.
            
            clear centrePoint
            clear fluenceWeightedDetector
            
            Approx_Outer_Angle = O * mradPerPixel % Displays outer angle into the command window.
            
            if size(meanSensitivity,2) < round(O*1.5) % Adds compatability for detectors that do not have their outer angle in the field of view.
                dataEnd = size(meanSensitivity,2) ;
            else
                dataEnd = round(O*1.5) ;
            end
            
            
            radialSensitivity = meanSensitivity(1:dataEnd)/mean(meanSensitivity(I:O)) ;
            
            flatnessScore =  100 * std(meanSensitivity(I:O)/mean(meanSensitivity(I:O))) ;
            
            clear I
            
            figure('Position',[108 148 fullscreen(3)-250 fullscreen(4)-396],'Name',horzcat('Detector Radial Sensitivity - ',version_string),'NumberTitle','off' ) ;
            bar(detectorAngleAxis(1:dataEnd),radialSensitivity,'FaceColor',[.361 0.749 0.749])
            xlim([0 Approx_Outer_Angle * 1.5])
            set(gca,'FontSize',28,'LineWidth',2,'XColor','k','Color',[0.914 0.914 0.914])
            set(gca,'layer','top')
            hline = refline([0 1]);
            set(hline,'Color','r','LineWidth',3)
            xlabel('Detector Angle','FontSize',34,'FontWeight','bold')
            ylabel('Relative Sensitivity','FontSize',34,'FontWeight','bold')
            
            % Detector angles expressed as a fraction of the inner angle:
            %                         figure('Position',[108 148 fullscreen(3)-250 fullscreen(4)-396],'Name',horzcat('Detector Radial Sensitivity - ',version_string),'NumberTitle','off' ) ;
            %                         bar(detectorAngleAxis(1:dataEnd)./innerAngle,radialSensitivity,'FaceColor',[1 0.616 0.435])
            %                         xlim([0 Approx_Outer_Angle/innerAngle * 1.5])
            %                         set(gca,'FontSize',28,'LineWidth',2,'XColor','k','Color',[0.914 0.914 0.914])
            %                         set(gca,'layer','top')
            %                         hline = refline([0 1]);
            %                         set(hline,'Color','r','LineWidth',3)
            %                         xlabel('Fraction of Inner Angle','FontSize',34,'FontWeight','bold')
            %                         ylabel('Relative Sensitivity','FontSize',34,'FontWeight','bold')
            
            clear dataEnd
            clear detectorAngleAxis
            clear AX
            clear O
            clear meanSensitivity
            clear deducedElectronDistribution
            clear weightedSensitivity
            
        else
            
            % Detector image is reshaped for histogram analysis:
            detectorHistogram = sort( reshape(detectorImage,1,[]) ) ;
            
            % View detector histogram:
            figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Detector Segmentation & Calibration - ',version_string),'NumberTitle','off' ) ;
            
            subplot(2,1,1)
            imagesc(detectorImage)
            axis image
            set(gca,'XTick',[],'YTick',[])
            title('Scanned Detector Map','FontSize',11,'FontWeight','bold')
            
            subplot(2,1,2)
            [N,binCenters] = hist(detectorHistogram,256) ;
            bar(binCenters,N);
            %         hold on
            %         plot(binCenters,gradient(N),'r')
            %         hold off
            xlim([min(min(detectorImage)) max(max(detectorImage))])
            grid on
            title('Detector Map Pixel Intensity Histogram','FontSize',11,'FontWeight','bold')
            xlabel('Detector Map Pixel-Intensities','FontWeight','bold')
            ylabel('Pixel Counts','FontWeight','bold')
            [a,~] = getpts(gcf) ;
            
            clear detectorHistogram
            
            vacuumIntesity = mean(detectorImage(detectorImage < a(1)))
            detectorIntesity = mean(detectorImage(detectorImage > a(2)))
            
            subplot(2,1,2)
            hBar = bar(binCenters,N);
            xlim([min(min(detectorImage)) max(max(detectorImage))])
            v = axis ;
            title('Detector Map Pixel Intensity Histogram','FontSize',11,'FontWeight','bold')
            xlabel('Detector Map Pixel-Intensities','FontWeight','bold')
            ylabel('Pixel Counts','FontWeight','bold')
            hold on
            line([a(1) a(1)],[0 v(4)],'LineWidth',3,'Color',[0 0 0]) ;
            line([a(2) a(2)],[0 v(4)],'LineWidth',3,'Color',[0 0 0]) ;
            line([vacuumIntesity vacuumIntesity],[0 v(4)],'LineStyle','--','LineWidth',2,'Color',[0 0 0]) ;
            line([detectorIntesity detectorIntesity],[0 v(4)],'LineStyle','--','LineWidth',2,'Color',[0 0 0]) ;
            hold off
            
            index = linspace( min(min(detectorImage)), max(max(detectorImage)) ,size(N,2)*5+1);
            coloring = ones(1,size(index,2)) ;
            coloring(index<a(1)) = 0 ;
            coloring(index>a(2)) = 2 ;
            
            bar_child=get(hBar,'Children');
            set(bar_child, 'CData',coloring);
            
            clear bar_child
            
            mycolor = [1 0 0;0 1 0;0 0 1]; % Specifies the colormap for the three segments of color from low to high.
            colormap(mycolor);
            clear mycolor
            
            %*****
            % Section generates three layer false colour image of
            % thresholded regions:
            
            thresholded_detector_image(:,:,1) = detectorImage < a(1) ;
            thresholded_detector_image(:,:,2) = 1 - (detectorImage > a(1)) - (detectorImage < a(2)) ;
            thresholded_detector_image(:,:,3) = detectorImage > a(2) ;
            thresholded_detector_image = single(thresholded_detector_image) ;
            
            subplot(2,1,1)
            imshow(cat(3,thresholded_detector_image(:,:,1),thresholded_detector_image(:,:,2),thresholded_detector_image(:,:,3)),'Border','tight')
            title('Scanned Detector Map','FontSize',11,'FontWeight','bold')
            axis on
            set(gca,'XTick',[],'YTick',[])
            
            clear thresholded_detector_image
            %*****
            
            beamIntensity = detectorIntesity - vacuumIntesity ;
            
            clear a
            clear index
            clear detectorImage
            %         clear detectorIntesity
            %         clear vacuumIntesity
        end
    end
    
    %% User loads STEM image:
    run 'C:\Users\Lewys\Box Sync\DPhil\MatLab Software\File Opener\File_Opener_1_5_5.m'
    if size(input_file,3) > 1
        warning('Line added at 428 to strip multiple frames. Analysing first frame only.')
%         input_file(:,:,1) = [] ;
        input_file(:,:,2:end) = [] ;
    end
    if loadSimulation == 1
        input_file = repmat(input_file,simulationTiling,simulationTiling) ;
    end
    
    input_file = imresize(input_file, upsamplingFactor) ;
    
    if addNoiseTargetSNR > 0
        minLog = min(min(input_file)) ;
        variance = ( std(reshape(input_file,1,[])) / addNoiseTargetSNR ) ^ 2 ;
        input_file = imnoise(input_file,'gaussian',0,variance ) ;
        clear variance
    end
    
    if loadSimulation == 0
        input_file = input_file - vacuumIntesity ;
        input_file = input_file / ( currentNormalisation * beamIntensity ) ; % The scaled ADF image. Units fraction of incident beam.
        
        clear currentNormalisation
        % This completes the intensity normalisation.
    end
    
    % Next the pixel width in units of nm should be determined from the image FT:
    if manualPixelWidth == 0
        figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Diffractogram Calibration - ',version_string),'NumberTitle','off' ) ;
        subplot(1,2,1)
        imagesc(input_file)
        set(gca,'XTick',[],'YTick',[])
        axis image
        colormap gray
        title('STEM Image','FontSize',11,'FontWeight','bold')
        
        % Generate the image FT:
        imageFT = log(medfilt2(abs(fftshift(fft2(input_file))))) ;
        imageFT(imageFT<0) = 0 ;
        
        subplot(1,2,2)
        imagesc(imageFT)
        set(gca,'XTick',[],'YTick',[])
        daspect([fliplr(size(imageFT)) 1]) % This line corrects the display aspect ratio of the Fourier transform.
        colormap gray
        title('Image Fourier Transform','FontSize',11,'FontWeight','bold')
        zoom(3*upsamplingFactor)
        
        [a,b] = getpts(gcf) ;
        center_point = [round((size(input_file,2)+1)/2) round((size(input_file,1)+1)/2)] ; % Coordinates of the center FT spot, notation [x y].
        
        pause(0.5)
        
        subplot(1,2,2)
        imagesc(imageFT)
        set(gca,'XTick',[],'YTick',[])
        daspect([fliplr(size(imageFT)) 1]) % This line corrects the display aspect ratio of the Fourier transform.
        colormap gray
        title('Image Fourier Transform','FontSize',11,'FontWeight','bold')
        zoom(3*upsamplingFactor)
        hold on
        line( [center_point(1) a(1)] , [center_point(2) b(1)] , 'Color', 'r', 'LineWidth' , 3 )
        line( [center_point(1) a(2)] , [center_point(2) b(2)] , 'Color', 'r', 'LineWidth' , 3 )
        hold off
        
        bigN = sum(calibrationPlanes.^2) ;
        planeSpacing = latticeParameter / sqrt(bigN) ; % The real-space separation of the planes in *nanometers*.
        
        vector_a     = [ b(1)-center_point(2) a(1)-center_point(1) ] ; % 1st basis vector.
        vector_b     = [ b(2)-center_point(2) a(2)-center_point(1) ] ; % 2nd basis vector.
        
        % Before solving these we take the absolute values as only the
        % magnitude of the pixel sizes are of interest:
        vector_a = abs(vector_a.^2) ;
        vector_b = abs(vector_b.^2) ;
        
        % Next these two vectors should be concatonated into a single array:
        vectors = vertcat(vector_a,vector_b) ;
        
        % And similarly for the vectors lengths:
        %         lengthVector_a = sqrt(sum(vector_a.^2)) ;
        %         lengthVector_b = sqrt(sum(vector_b.^2)) ;
        %         lengths = vertcat(lengthVector_a,lengthVector_b) ;
        
        % Then the system of equations can be solved. Units of reciprocal
        % nanometers.
        reciprocal_pixel_dimensions = realsqrt( abs ( vectors \ [ 1/planeSpacing ; 1/planeSpacing ].^2 ) ) ; % Has units of "the number of reciprocal nm per FT pixel".
        
        % Which can then be expressed in real-space nanometers.
        real_pixel_dimensions = 1 ./ reciprocal_pixel_dimensions ./ size(imageFT)' ;
        
        clear imageFT
        clear center_point
        clear vector_a
        clear vector_b
        
        pixelWidth = real_pixel_dimensions(1,1) ; % pixelWidth here has units of nanometers.
        
        pixelArea = real_pixel_dimensions(1,1) * real_pixel_dimensions(2,1) ; % pixelArea here has units of square nanometers.
        
        % Variable clean up:
        clear a
        clear b
        clear lengthVector_a
        clear lengthVector_b
        clear bigN
        clear calibrationPlanes
        clear planeSpacing
        clear latticeParameter
        clear averageVectorLength
        
        % Generate scale bar for real-space image:
        targetLength = (0.5 * round(( round ( size(input_file,1) / 5 ) * pixelWidth ) * 2 )) ; % Creates a scalebar roughly one fifthe the width of the image rounded to one decimal place of nanometres.
        scalebarLengthPixels =  targetLength / pixelWidth ;
        startCorner = round(size(input_file,1) / 10) ;
        
        pause(0.5)
        
        subplot(1,2,1)
        imagesc(input_file)
        set(gca,'XTick',[],'YTick',[])
        axis image
        colormap gray
        title('STEM Image','FontSize',11,'FontWeight','bold')
        hold on
        line( [startCorner startCorner+scalebarLengthPixels] , [startCorner startCorner] , 'Color', 'r', 'LineWidth' , 3 )
        text(startCorner+0.3*scalebarLengthPixels , startCorner-0.1*scalebarLengthPixels, horzcat(num2str(targetLength),'nm'),'FontSize',11,'FontWeight','bold','Color','r')
        hold off
        
        pause(0.5)
        
        clear targetLength
        clear scalebarLengthPixels
        clear startCorner
    else
        pixelWidth = manualPixelWidth / upsamplingFactor ;
        pixelArea = pixelWidth^2 ;
    end
    
    disp('Loading pre-determined peak-poitions.')
    %     load('C:\Users\Lewys\Box Sync\PostDoc\ESTEEM2\Noise Study\December 2014 Tests\Distorted Data Sets\Distortion Free\Perfect_points.mat')
%     load('C:\Users\Lewys\Box Sync\PostDoc\ESTEEM2\Noise Study\December 2014 Tests\Rotated Distorted Data Sets\Perfect_peaks.mat')
%     total_atoms = size(point_coordinates,1) ;
    
    % END preloading lines.
    
    if exist('point_coordinates' , 'var' ) == 0
        run 'C:\Users\Lewys\Box Sync\DPhil\MatLab Software\Ranger\Ranger_Core_2_4.m'
    end
    
    if addRandomOffset > 0
        offsets =  addRandomOffset * 2 * ( rand(size(point_coordinates)) - 0.5 ) ;
        point_coordinates = point_coordinates + offsets ;
        StdOffsets = mean(std(offsets)) ;
        clear offsets
        
        point_coordinates_integer = round(point_coordinates) ;
        for i = 1:size(point_coordinates_integer,1)
            peak_values(i) = input_file(point_coordinates_integer(i,2),point_coordinates_integer(i,1)) ;
        end
        
        StdPeaks = std(peak_values) / mean(peak_values) ;
        clear peak_values
        clear StdPeaks
        
    end
    
    distance_log = zeros(total_atoms,1) ;
    % Next the distance to features must be calcualted for each pixel:
    
    % Generates the measurement progress bar:
    h2=waitbar(0,{'Now Determining Voroinoi Cells.'},'Name',version_string);
    hw2=findobj(h2,'Type','Patch');
    set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % changes the color to red.
    
    if autoRadius == 1
        integrationRadius = best_size * 1.2 ; % Modifier to limit total inclusion radius.
    else
        integrationRadius = manualRadius ;
    end
    
    clear manualRadius
    
    point_record = zeros(size(input_file)) ;
    
    for i = 1:size(search_record,1)
        for j = 1:size(search_record,2)
            % For every pixel the distance to all points must be calcualted:
            
            distance_log = sqrt( (point_coordinates(:,2)-i).^2 + (point_coordinates(:,1)-j).^2 ) ;
            
            % Next for that pixel the minimum distance to and point should be
            % checked and discarded if too large:
            [distMin, minIndex] = min(distance_log) ;
            if distMin > integrationRadius
                point_record(i,j) = 0 ;
            else
                point_record(i,j) = minIndex ;
            end
        end
        
        if (i/10) == round(i/10)
            percentageDetermined = i ./ size(search_record,1) ;
            waitbar_progress_text1 = horzcat ( 'Now Determining Voroinoi Cells.' ) ;
            waitbar_progress_text2 = horzcat ( 'Progress = ' , num2str(100*percentageDetermined) , '%.' ) ;
            waitbar(percentageDetermined,h2,{waitbar_progress_text1;waitbar_progress_text2},'Name',version_string) ;
        end
    end
    
    close (h2)
    clear i
    clear j
    clear point
    clear percentageDetermined
    
    % figure
    % imagesc(point_record)
    % axis image
    % colormap (prism(total_atoms))
    % title ('Identified Feature Areas')
    % hold on
    % plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
    %     'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
    % hold off
    
    % Next for each of these a mask can be defined and the integrated intensity
    % output:
    
    integratedIntensity = zeros(total_atoms,1) ;
    if subtractBackground == 1
        backgroundIntensity = zeros(total_atoms,1) ;
    end
    featureArea = zeros(total_atoms,1) ;
    
    % Up to this point the greyscale of the image as as a fraction of incident
    % probe. The DC offset has already been subtracted but it may also be
    % necessary to subtract a background from carbon for example.
    
    if subtractBackground == 1
        if fitLocalBackground == 1
            
            
            particleMask = ceil(search_record) ;
            se = strel('disk',round(best_size*1.5),0);
            particleMask = imdilate(particleMask,se);
            clear se
            
            if add_extra_mask == 1
                title('Please draw on extra mask. 1) Left-click to define corners, double-click to close loop. 2) Left-click to drag corners/shape, double-click to accept.','FontSize',11,'FontWeight','bold')
                extra_mask = roipoly ;
                particleMask = particleMask + extra_mask ;
            end
            
            particleMask = particleMask / max(particleMask(:)) ;
            particleMask = ceil(particleMask) ;
            
            
            backgroundImage = roifill(input_file, particleMask) ; % We can use this as an initial estimate of the background image.
            
            startingBackgroundGuess = backgroundImage ;
            
            sigma = round(0.33*best_size) ;    % The standard deviation (width) of the image smoothing kernel.
            x = -3*sigma:3*sigma ; % Create canvas for the 1D Gaussian.
            G = 1 / (sigma * sqrt(2 * pi)) * exp(-x.^2 / (2*sigma^2)); % 1D Gaussian
            G2 = G' *  G; % 2D
            G2 = G2 / sum(G2(:)) ;
            clear G
            
            backfill_fig = figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Local Background Subtraction - ',version_string),'NumberTitle','off' ) ;
            subplot(1,2,1)
            imagesc(startingBackgroundGuess)
            axis image
            set(gca,'xtick',[],'ytick',[])
            title('Current Background Estimate', 'FontWeight','Bold','FontSize',11 )
            v = caxis ; % Measures the colour scaling of the image with particle.
            
            % Generates the measurement progress bar:
            h2=waitbar(0,'Attempting Iterative Background Recreation','Name',version_string,'CreateCancelBtn','abort = 1 ;');
            abort = 0 ;
            hw2=findobj(h2,'Type','Patch');
            set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % changes the color to red.
            
            convergenceMetric = zeros(1,backfillIterations) ;
            for backfillIteration = 1:backfillIterations
                
                currentBackgroundGuess = padarray(startingBackgroundGuess, [ 3*sigma 3*sigma ] ,'replicate','both') ;
                %             currentBackgroundGuess = medfilt2(currentBackgroundGuess) ;
                %             currentBackgroundGuess = conv2(currentBackgroundGuess,G2,'valid') ;
                currentBackgroundGuess = convolve2(currentBackgroundGuess,G2,'valid') ; % 2D convolution implimented usign the University Sussex method Copyright David Young. Alternative is the inbuilt 'conv2'.
                
                rebuiltBackground = (particleMask==0).*input_file + (particleMask>0).*currentBackgroundGuess ;
                
                subplot(1,2,1)
                imagesc(rebuiltBackground)
                axis image
                set(gca,'xtick',[],'ytick',[])
                title('Current Background Estimate', 'FontWeight','Bold','FontSize',11 )
                caxis(v)
                
                convergenceMetric(backfillIteration) = sum(sum(abs(rebuiltBackground - startingBackgroundGuess))) / sum(sum(particleMask>0))  ;
                
                subplot(1,2,2)
                semilogy(convergenceMetric/max(convergenceMetric(:)))
                title('Relative Change in Background', 'FontWeight','Bold','FontSize',11 )
                xlabel('Iteration Number', 'FontWeight','Bold','FontSize',11 )
                ylabel('Iteration Residual', 'FontWeight','Bold','FontSize',11 )
                axis square
                pause(0.01)
                
                startingBackgroundGuess = rebuiltBackground ;
                
                waitbar(backfillIteration/backfillIterations,h2,'Attempting Iterative Background Recreation','Name',version_string) ;
                
                if abort == 1
                    clear abort
                    break ;
                end
                
            end
            
            delete(h2)
            close(backfill_fig)
            
            backgroundImage = rebuiltBackground ; % After the iterattive refinement the background image is updated to the last loop's output.
            
            clear startingBackgroundGuess
            clear backfillIterrations
            clear rebuiltBackground
            clear currentBackgroundGuess
            clear particleMask
            clear G2
            
            figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Local Background Subtraction - ',version_string),'NumberTitle','off' ) ;
            
            subplot(2,1,1)
            imagesc(input_file)
            title('Experimental Image', 'FontWeight','Bold','FontSize',11 )
            set(gca,'xtick',[],'ytick',[])
            axis image
            colorbar
            v = caxis ; % Reasures the colour scaling of the image with particle.
            
            subplot(2,2,3)
            imagesc(backgroundImage)
            title('Background Image', 'FontWeight','Bold','FontSize',11 )
            set(gca,'xtick',[],'ytick',[])
            axis image
            caxis(v)
            colorbar
            
            subplot(2,2,4)
            imagesc(input_file - backgroundImage)
            title('Sample Image', 'FontWeight','Bold','FontSize',11 )
            set(gca,'xtick',[],'ytick',[])
            axis image
            colorbar
            
            clear v
            
            input_file = input_file - backgroundImage ;
            
            
            %% Developemt analysis figures - all code redundant:
            % imageHistogram = reshape(input_file,1,[]) ;
            %
            % figure
            % subplot(2,1,1)
            % hist(imageHistogram,256)
            % axis square
            %
            % particleMask = ceil(search_record) ;
            % se = strel('disk',best_size,0);
            % particleMask = imdilate(particleMask,se);
            % clear se
            %
            % imageHistogram = nonzeros(particleMask .* backgroundImage) ;
            % reshape(backgroundImage,1,[]) ;
            % subplot(2,2,3)
            % hist(imageHistogram,256)
            % axis square
            %
            % particleMask = ceil(search_record) ;
            % se = strel('disk',best_size,0);
            % particleMask = imdilate(particleMask,se);
            % clear se
            %
            % imageHistogram = nonzeros(input_file .* particleMask ) ;
            % subplot(2,2,4)
            % hist(imageHistogram,256)
            % axis square
            % clear imageHistogram
            
            clear particleMask
            %             clear backgroundImage
            
            
        else
            Hist = figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Global Background Selection Dialogue - ',version_string),'NumberTitle','off' ) ;
            hist(reshape(input_file,1,[]), 256 ) % Display distribution of peak intensities.
            title( {'Intensity Profile of Image Pixels';'Click centre of background peak.';'"Backspace" to undo, "Return" to accept.'} , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
            xlabel('Feature Intensity')
            ylabel('Number of Features')
            [meanBackground,~] = getpts(gcf) ;
            close(Hist)
            input_file = input_file - meanBackground ;
            pause(0.1)
        end
    end
    
    % Generates the measurement progress bar:
    h2=waitbar(0,{'Now Integrating Voroinoi Cells.'},'Name',version_string);
    hw2=findobj(h2,'Type','Patch');
    set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % changes the color to red.
    
    if calcualte_downsampled_XS == 1
        % First we generate a diminished version of "point_record"
        [downsampled_X,downsampled_Y] = meshgrid(1:size(point_record,1), 1:size(point_record,2)) ;
        downsampled_X = mod(downsampled_X,2) ;
        %         downsampled_X = 1 - downsampled_X ; % Debug line only: inverts checkerboard for downsampling.
        downsampled_Y = mod(downsampled_Y,2) ;
        %         downsampled_Y = 1 - downsampled_Y ; % Debug line only: inverts checkerboard for downsampling.
        downsampling_template = downsampled_X .* downsampled_Y ;
        point_record_downsampled = downsampling_template .* point_record ;
        clear downsampled_X
        clear downsampled_Y
        clear downsampling_template
        integratedIntensity_downsampled = zeros(total_atoms,1) ;
        featureArea_downsampled = zeros(total_atoms,1) ;
    end
    
    for point = 1:total_atoms
        currentMask = abs (point_record == point) ;
        currentFeature = currentMask .* input_file ;
        integratedIntensity(point) = sum(currentFeature(currentFeature ~= 0)) ;
        if subtractBackground == 1 && fitLocalBackground == 1
            currentFeature = currentMask .* backgroundImage ;
            %             backgroundVaraincePerCell(point) = (std(nonzeros(currentFeature))).^2 ;
            backgroundIntensity(point) = sum(currentFeature(currentFeature ~= 0)) ;
        end
        featureArea(point) = nnz(currentFeature) ;
        
        if calcualte_downsampled_XS == 1
            currentMask = abs (point_record_downsampled == point) ;
            currentFeature = currentMask .* input_file ;
            integratedIntensity_downsampled(point) = sum(currentFeature(currentFeature ~= 0)) ;
            featureArea_downsampled(point) = 4 * nnz(currentFeature) ; % Note the '4x' factor here is for the assumption of a larger pixel size on downsampling.
        end
        
        if mod(point,10) == 0
            waitbar_progress_text1 = horzcat ( 'Now Integrating Voroinoi Cells.' ) ;
            waitbar_progress_text2 = horzcat ( 'Progress = ' , num2str(100*(point/total_atoms)) , '%.' ) ;
            waitbar((point/total_atoms),h2,{waitbar_progress_text1;waitbar_progress_text2},'Name',version_string) ;
        end
    end
    
    delete(h2)
    clear currentMask
    clear currentFeature
    
    % These integrted intensities are then scaled to give cross-sections:
    pixelScaling = ( pixelArea ) * 10.^4 ; % Units Mb
    integratedIntensity = integratedIntensity * pixelScaling ; % Conversion to cross-sections.
    if calcualte_downsampled_XS == 1
        integratedIntensity_downsampled = integratedIntensity_downsampled * 4 * pixelScaling ; % Conversion to cross-sections. Note '4x' factor to allow for downsampling.
    end
    
    % Once the integrated intensities are know these can be rebuilt into an
    % image corresponding to the original positions:
    intensityRecord = zeros(size(input_file)) ;
    for i = 1 : size(point_record,1)
        for j = 1 : size(point_record,2)
            if point_record(i,j) > 0
                intensityRecord(i,j) = integratedIntensity(point_record(i,j)) ;
            else
                intensityRecord(i,j) = 0 ;
            end
        end
    end
    
    figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Cross-section Map - ',version_string),'NumberTitle','off' ) ;
    imagesc(intensityRecord)
    colorbar
    axis image
    set(gca,'XTick',[],'YTick',[])
    title ('Identified Feature Cross Sections (Mb)' , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
    hold on
    plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
        'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
    hold off
    
    % Next remove zeros:
    disp('Warning: experiemnatal zero removal suppressed!')
%     featureArea = featureArea(featureArea ~= 0) ;
%     integratedIntensity = integratedIntensity(featureArea ~= 0) ;
%     if calcualte_downsampled_XS == 1
%         featureArea_downsampled = featureArea_downsampled(featureArea ~= 0) ;
%         integratedIntensity_downsampled = integratedIntensity_downsampled(featureArea ~= 0) ;
%     end
    
    figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Integration Zone Area Histogram - ',version_string),'NumberTitle','off' ) ;
    hist(featureArea,total_atoms/5)
    xlabel('Integration Area')
    ylabel('Histogram Counts')
    if calcualte_downsampled_XS == 1
        set(get(gca,'child'),'FaceColor','none','EdgeColor','r');
        hold on
        hist(featureArea_downsampled,total_atoms/10) ;
        hold off
    end
    
    % maxArea = getpts(gcf) ;
    % integratedIntensity(featureArea>maxArea) = [] ;
    
    figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Scattering Cross-section Histogram - ',version_string),'NumberTitle','off' ) ;
    %     [nelements,xcenters] = hist(integratedIntensity,total_atoms/16) ;
    hist(integratedIntensity,total_atoms/10) ;
    xlabel('Scattering Cross Sections (Mb)', 'FontWeight' , 'Bold' , 'FontSize' , 11)
    ylabel('Histogram Counts', 'FontWeight' , 'Bold' , 'FontSize' , 11)
    if calcualte_downsampled_XS == 1
        set(get(gca,'child'),'FaceColor','r','EdgeColor','k');
        hold on
        hist(integratedIntensity_downsampled,total_atoms/10) ;
        hold off
        legend('Integrated Cross-sections','Downsampled (2x2) Cross-sections')
    end
    % xlim([0 1.5e8])
    
    if calcualte_downsampled_XS == 1
        figure
        plot(integratedIntensity,integratedIntensity_downsampled,'b.')
        title('Cross-section Variance Analysis' , 'FontWeight' , 'Bold' , 'FontSize' , 13 )
        xlabel('Integrated Cross-secction (Mb)' , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
        ylabel('Downsampled (2x2) Cross-section (Mb)' , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
        axis square
        grid on
    end
    % StdPeaks = std(peak) / mean(peak)
    
    if subtractBackground == 1
        figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Scattering Cross-section Histogram - ',version_string),'NumberTitle','off' ) ;
        hist(backgroundIntensity,total_atoms/10) ;
        xlabel('Background Intensity')
        ylabel('Histogram Counts')
    end
    
    MeanCrossSection = mean(integratedIntensity)
    StdCrossSection = std(integratedIntensity)
    
    if addNoiseTargetSNR > 0
        contrast = mean(peak) / minLog
    else
        
    end
    
    
    
    %% Experimental section for Gaussian asymmetry
    
    if calculateEllipticity == 1
        
        figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Cross-section versus Ellipticity Plot - ',version_string),'NumberTitle','off' ) ;
        
        x = 0 ;
        h = waitbar(x,'Calculating Column Ellipticities') ;
        
        height_A = zeros( 1 , size(point_coordinates,1) ) ;
        sigma_a = zeros( 1 , size(point_coordinates,1) ) ;
        sigma_b = zeros( 1 , size(point_coordinates,1) ) ;
        sigma_c = zeros( 1 , size(point_coordinates,1) ) ;
        offset_o = zeros( 1 , size(point_coordinates,1) ) ;
        ramp_r1 = zeros( 1 , size(point_coordinates,1) ) ;
        ramp_r2 = zeros( 1 , size(point_coordinates,1) ) ;
        centre_x0 = zeros( 1 , size(point_coordinates,1) ) ;
        centre_y0 = zeros( 1 , size(point_coordinates,1) ) ;
        rebuilt_image = zeros(size(input_file)) ;
        
        [Yall,Xall]=meshgrid(1:size(input_file,2),1:size(input_file,1)) ;   % Generates a mesh with those dimensions.
        
        for feature = 1 : size(point_coordinates,1)
            
            X = Xall - round(point_coordinates(feature,2)) ;
            Y = Yall - round(point_coordinates(feature,1)) ;
            
            % Calculate a circular mask to put around the each feature:
            mask = 1*(point_record == feature) ;
            
            Z = input_file.*mask ;
            
            [fitresult, gof] = weightedEllipticalGaussianFit(Y, X, Z, mask) ;
            
            height_A(feature) = fitresult.A ;
            sigma_a(feature)  = fitresult.a ;
            sigma_b(feature)  = fitresult.b ;
            sigma_c(feature)  = fitresult.c ;
            offset_o(feature) = fitresult.o ;
            ramp_r1(feature)  = fitresult.r1 ;
            ramp_r2(feature)  = fitresult.r2 ;
            centre_x0(feature) = fitresult.x0 ;
            centre_y0(feature) = fitresult.y0 ;
            
            rebuilt_image = rebuilt_image + mask.*(fitresult.A*exp( - (fitresult.a*(Y-fitresult.y0).^2 + 2*fitresult.b*(Y-fitresult.y0).*(X-fitresult.x0) + fitresult.c*(X-fitresult.x0).^2)) + fitresult.r1*Y + fitresult.r2*X + fitresult.o )  ;
            
            imagesc(rebuilt_image)
            axis image
            colorbar
            
            pause(0.05)
            
            x = feature / size(point_coordinates,1) ;
            waitbar(x)
            
        end
        close(h)
        clear mask
        
        sigma_x = abs(sqrt( sigma_a + sigma_c + sqrt((4*sigma_b.^2)+(sigma_a-sigma_c).^2) )) ;
        sigma_y = abs(sqrt( sigma_a + sigma_c - sqrt((4*sigma_b.^2)+(sigma_a-sigma_c).^2) )) ;
        theta   = 0.5 * atand( (2*sigma_b) ./ (sigma_c - sigma_a) ) ;
        
        ellipticity = (sigma_x ./ sigma_y) - 1 ;
        
        clear height_A
        clear sigma_a
        clear sigma_b
        clear sigma_c
        clear offset_o
        clear ramp_r1
        clear ramp_r2
        clear centre_x0
        clear centre_y0
        
        clear sigma_x
        clear sigma_y
        
        meanEllipticity = mean(ellipticity)
        stDevEllipticity = std(ellipticity)
        
        figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Cross-section versus Ellipticity Plot - ',version_string),'NumberTitle','off' ) ;
        plot(integratedIntensity,ellipticity,'o')
        xlabel('Scattering Cross-section (Mb)')
        ylabel('Ellipticity')
        
        quiver_x = sind(theta)'.*ellipticity' ;
        quiver_y = -cosd(theta)'.*ellipticity' ;
        
        figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat('Ellipticity Map - ',version_string),'NumberTitle','off' ) ;
        imagesc(rebuilt_image)
        axis image
        colorbar
        set(gca,'XTick',[],'YTick',[])
        title ('Estimated Asymmetries' , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
        hold on
        quiver(point_coordinates(:,1),point_coordinates(:,2),quiver_x,quiver_y,'-k','LineWidth',2)
        hold off
        
    end
    
    if use_SpectrumIntegrator == 1
        % User loads spectral data:
        h = msgbox('Loading spectral data. May take > 1min for large files.',version_string,'help') ;
        run 'C:\Users\Lewys\Box Sync\DPhil\MatLab Software\File Opener\File_Opener_1_5_5'
        close(h)
        spectral_data = input_file ;
        clear input_file
        
        integrated_spectrum = zeros(total_atoms , size(spectral_data,3) ) ;
        for point = 1:total_atoms
            currentMask = abs (point_record == point) ;
            for dispersion = 1:size(spectral_data,3)
                currentFeature = currentMask .* spectral_data(:,:,dispersion) ;
                integrated_spectrum(point,dispersion) = sum(sum(currentFeature(currentFeature ~= 0))) ;
            end
            
            % Put some kind of waitbar here!!
            
        end
        
        % SpectrumIntegrator display:
        figure
        subplot(2,1,1)
        imagesc(intensityRecord)
        colorbar
        axis image
        set(gca,'XTick',[],'YTick',[])
        title ('Identified Feature Cross Sections (Mb)' , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
        hold on
        plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
            'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
        hold off
        
        % Get location of cursor within loop:
        exit = 0 ;
        while exit == 0
            % Create draggable cursor:
            h = impoint(gca,[]);
            
            pos = round(getPosition(h)) ;
            feature_number = point_record(pos(2), pos(1)) ;
            
            if feature_number == 0
                disp('Click is out of bounds');
            else
                subplot(2,1,2)
                %             plot( [1/910 1/720 1/580 1/460 1/360 1/230 1/185 1/145 1/115 1/91 1/73 1/58 1/46 1/37 1/29.5] , integrated_spectrum(feature_number,:))
                plot(  integrated_spectrum(feature_number,:))
                title ('Spectral Profile' , 'FontWeight' , 'Bold' , 'FontSize' , 11 )
            end
            
            exit = waitforbuttonpress ;
            delete(h)
        end
        
        clear exit
        
    end
    
    
end

