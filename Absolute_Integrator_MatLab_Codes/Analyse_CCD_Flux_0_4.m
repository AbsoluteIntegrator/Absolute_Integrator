
% Load the CCD flux pattern:
run 'C:\Users\Lewys\Box Sync\DPhil\MatLab Software\File Opener\File_Opener_1_5_5'
flux_pattern = input_file ;
clear input_file

% Remove any negative pixel values:
flux_pattern(flux_pattern<0) = 0 ;
scaled_flux = flux_pattern / max(max(medfilt2(flux_pattern))) ;
clear flux_pattern

% Find the centre of the image:
mask = round(scaled_flux+0.25) ; % Creates a mask for all areas above 25% of image max.
s = regionprops(mask, 'centroid');
[cx] = round(s.Centroid(1)) ;
[cy] = round(s.Centroid(2)) ;

% [cy,cx] = find(mask, 1,'first') ;
% warning('Hard code here')
% cx = 1025
% cy = 973
clear mask

% Find distane to nearest edge:
distance_to_edges = [ cx , cy , size(scaled_flux,1)-cy , size(scaled_flux,2)-cx ] ;
distance_to_nearest_edge = min(distance_to_edges) ;
clear distance_to_edges

% Display the pattern with linear and log scales:
figure('Position',[108 148 fullscreen(3)-250 fullscreen(4)-396],'Name',horzcat('Detector Azimuthal Sensitivity - ',version_string),'NumberTitle','off' ) ;
pause(0.05) % Allow figure to update.

subplot(2,2,1)
imagesc(medfilt2(scaled_flux))
axis image
set(gca,'XTick',[],'YTick',[])
title('CCD Flux Pattern (linear scale)','FontSize',11,'FontWeight','bold')
hold on
plot(cx,cy,'p','MarkerEdgeColor','w','MarkerFaceColor','w') ;
hold off

subplot(2,2,2)
imagesc(log(medfilt2(scaled_flux)))
axis image
set(gca,'XTick',[],'YTick',[])
title('CCD Flux Pattern (log scale)','FontSize',11,'FontWeight','bold')
hold on
plot(cx,cy,'p','MarkerEdgeColor','w','MarkerFaceColor','w') ;
hold off

pause(0.05) % Allow figure to update.

% Calculate flux in new coordinate set:
linear_flux = zeros(360,distance_to_nearest_edge) ;
i     = 1;

for r=0:distance_to_nearest_edge-1 % !! Note radius starts at zero here.
    j = 1;
    for a=0:2*pi/360:2*pi-2*pi/360
        linear_flux(j,i) = scaled_flux(cy+round(r*sin(a)),cx+round(r*cos(a))) ; % New flux distribution is extracted pixel by pixel using a nearest neighbour type look-up.
        j = j + 1;
    end
    i = i + 1;
end
clear scaled_flux
clear a
clear r
clear i
clear j

subplot(4,2,5)
imagesc(linear_flux)
set(gca,'XTick',[],'YTick',[])
title('CCD Flux Pattern (linear scale)','FontSize',11,'FontWeight','bold')

subplot(4,2,6)
imagesc(log(linear_flux))
set(gca,'XTick',[],'YTick',[])
title('CCD Flux Pattern (log scale)','FontSize',11,'FontWeight','bold')

[pix_flux_cutoff, ~] = getpts ;

pause(0.05) % Allow figure to update.

% Calcualate position of edge of BF disk:
flux_gradient = gradient(mean(linear_flux,1)) ;
mrad_per_pixel_flux = convergence_angle / find(flux_gradient == min(flux_gradient)) ; % returns the number of pixels equal to the edge of the BF disk.
clear flux_gradient

flux_angle_axis = mrad_per_pixel_flux*(0:size(linear_flux,2)-1) ;
flux_cuttof_angle = pix_flux_cutoff  * mrad_per_pixel_flux ;

subplot(4,2,7)
plot(flux_angle_axis,mean(linear_flux,1) , 'LineWidth',2) ;
hold on
line([innerAngle innerAngle] , [0 1] , 'LineStyle' , '--' , 'Color' , 'g' , 'LineWidth',2)
line([flux_cuttof_angle flux_cuttof_angle] , [0 1] , 'LineStyle' , '--' , 'Color' , 'r' , 'LineWidth',2)
hold off
grid on
xlim([1 max(flux_angle_axis)])
xlabel('Scattering Angle (mrad)','FontWeight','bold')
ylabel('Mean Flux (linear scale)','FontWeight','bold')

subplot(4,2,8)
semilogy(flux_angle_axis,mean(linear_flux,1) , 'LineWidth',2) ;
axis tight
yL = get(gca,'YLim');
hold on
line([innerAngle innerAngle] , yL , 'LineStyle' , '--' , 'Color' , 'g' , 'LineWidth',2)
line([flux_cuttof_angle flux_cuttof_angle] , yL , 'LineStyle' , '--' , 'Color' , 'r' , 'LineWidth',2)
hold off
grid on
xlabel('Scattering Angle (mrad)','FontWeight','bold')
ylabel('Mean Flux (log scale)','FontWeight','bold')

% END PLOTTING
clear distance_to_nearest_edge
clear cx
clear cy
clear yL

% Determine region of plot to fit to power law:
start_pix = round( innerAngle / mrad_per_pixel_flux ) ;
flux_over_detector = mean(linear_flux(:,start_pix:pix_flux_cutoff),1) ;

% Perform fitting:
[xData, yData] = prepareCurveData( flux_angle_axis(start_pix:pix_flux_cutoff) , flux_over_detector );

% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [100 -3 0];

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );
clear opts
clear ft
clear start_pix


% Plot fit with data.
hold on
semilogy(xData,fitresult.a*xData.^fitresult.b+fitresult.c , 'c-' , 'LineWidth',2) ;
hold off
legend( 'Experimental Flux Profile', 'Detector Inner Angle' , 'Flux Cut-off' , 'Fitted Flux Profile' , 'Location', 'NorthEast' );

% Return TDS exponent:
TDSexponent = -fitresult.b ;

clear flux_over_detector
clear linear_flux
clear flux_angle_axis
clear fitresult
clear xData
clear yData
