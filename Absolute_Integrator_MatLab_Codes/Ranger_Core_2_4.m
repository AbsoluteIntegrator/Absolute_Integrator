clear file_type

% This entire core file is configured as a testing copy and will not run
% beyond a pre-determined expiration date. The following condition checks
% for the system date and allows a finite use period only:

if exist('now','var') == 1  % This line prevens circumvention of date check by the user specifying a variable called 'now'.
    clear now
end

run Ranger_license_check_1_0.m

fullscreen = get(0,'ScreenSize') ;
inputUnfiltered = input_file ;

[m,n] = size(input_file) ;

% **NEW** removed hard-coding of big!!
% First check smallest dimension of image:
smallestDim = min([m n]) ;
% Big must be around 1/8th of this. Any larger has no physcial meaning
% (fewer than 4 periodic features at the Nyquist limit). Big must also
% be rounded to the neareast ODD number.
big = 2 .* floor((smallestDim/8)/2) - 1 ;

clear smallestDim

if exist('best_size','var') == 0
    
    if show_progress == 1
        progFig = figure('Position',[108 148 fullscreen(3)-215 fullscreen(4)-336],'Name',horzcat(version_string,' Analysis Progress'),'NumberTitle','off' ) ;
    end
    
    size_axis = (3:2:big) ; % Introduce and define the range of trial feature spacings to investigate.
    k = NaN*ones(1,round((big-3)/2) + 1) ;
    
    % This first loop probes the number of responses as a function of
    % test box size:
    
    h2=waitbar(0,'Determining Optimum Peak-Separation.','Name',version_string,'CreateCancelBtn','abort = 1 ;');
    abort = 0 ;
    hw2=findobj(h2,'Type','Patch');
    set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % Changes the color to red.
    
    for trialSize = 3:2:big
        
        
        
        %% BEGIN CODE FOR FEATURE IDENTIFICATION - SHOULD REMAIN IDENTICAL TO COPY BELOW.
        inputOffset = single(input_file - min(input_file(:))) ;         % Removes 'DC offset' from image to simplify Gaussian fitting.
        vert_offset = zeros(m,n) ;                              % Create blank arrays.
        horz_offset = zeros(m,n) ;
        peak        = zeros(m,n) ;
        spread      = zeros(m,n) ;
        test_box_padding = ( trialSize - 1 ) / 2 ;              % Half of the trial size, equivalent to the border that will not be inspected.
        baseAxis = ( -test_box_padding : test_box_padding )  ;  % Coordinate set for X and Y fitting.
        baseAxis = baseAxis(:) ;
        A = [baseAxis.^2 , baseAxis , ones(size(baseAxis))] ;
        % Followed by the restoration progress bar:
        % h2=waitbar(0,'Identifying Image Peaks...','Name',version_string);
        % hw2=findobj(h2,'Type','Patch');
        % set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % Changes the color to red.
        for i = test_box_padding + 1 : m - ( test_box_padding + 1 )
            
            currentStrip = inputOffset( i - test_box_padding : i + test_box_padding , : ) ;
            
            for j = test_box_padding + 1 : n - ( test_box_padding + 1 )
                I = currentStrip( : ,  j - test_box_padding : j + test_box_padding ) ;
                
                horizontalProfile = sum(I,1) ; % Generate two 1D profiles for Gaussian fitting.
                b = log(horizontalProfile(:)); % Natural log of the 1D data
                solutionH = A \ b ;                 % Least-squares solution
                horz_offset(i,j) = -solutionH(2)/solutionH(1)/2 ; % The horizontal offset refinement.
                
                verticalProfile = sum(I,2)' ; % Generate two 1D profiles for Gaussian fitting.
                b = log(verticalProfile(:)); % Natural log of the 1D data
                solutionV = A \ b ;                 % Least-squares solution
                vert_offset(i,j) = -solutionV(2)/solutionV(1)/2 ; % The vertical offset refinement.
                
                peak(i,j) = max(horizontalProfile(:)+verticalProfile(:))/(2*trialSize) ; % Calculates the height of the fitted Gaussian.
                spread(i,j) = real(sqrt( -1/2/(0.5*solutionH(1)+0.5*solutionV(1)) )) ; % Calculates the width of the fitted Gaussian.
            end
            
            percentageRefined = ( ((trialSize-3)/2) / ((big-1)/2) ) +   ( ( (i-test_box_padding) / (m - 2*test_box_padding) ) / (((big-1)/2))) ; % Progress metric when using a looping peak-finding waitbar.
            waitbar(percentageRefined,h2) ;
        end
        clear A
        clear b
        clear solutionH
        clear solutionV
        % delete (h2)
        % Feature identification section:
        normalisedPeak = (peak ./ ( max(max(inputOffset)) - min(min(inputOffset))) ) ; % Make use of peak knowledge:
        normalisedPeak(normalisedPeak<0) = 0 ; % Forbid negative (concave) Gaussians.
        spread = spread ./ trialSize ;         % Normalise fitted Gaussian widths.
        offsetRadius = realsqrt( (horz_offset).^2 + (vert_offset).^2 ) ; % Calculate offset radii.
        offsetRadius = offsetRadius ./ trialSize ;
        offsetRadius(offsetRadius==0) = 0.001 ; % Remove zeros values to prevent division error later.
        % Create search metric and screen impossible peaks:
        search_record = ( normalisedPeak ./ offsetRadius ) ;
        search_record(search_record > 1) = 1 ;
        search_record(search_record < 0) = 0 ;
        search_record(spread < 0.05) = 0 ;      % Invalidates negative Gaussian widths.
        search_record(spread > 1) = 0 ;         % Invalidates Gaussian widths greater than a feature spacing.
        search_record(offsetRadius > 1) = 0 ;   % Invalidates Gaussian widths greater than a feature spacing.
        search_record = medfilt2(search_record, [round(trialSize/3) round(trialSize/3)]) ; % Median filter to strip impossibly local false-positive features.
        sensitivityThreshold = 0.34 ;           % This is an Admin tunable parameter that is defined here within the core file.
        search_record(search_record < sensitivityThreshold ) = 0 ;  % Collapse improbable features to zero likelyhood.
        search_record(search_record >= sensitivityThreshold ) = 1 ; % Round likelyhood of genuine features to unity.
        search_record = bwmorph(search_record,'shrink',Inf) ;       % Errod regions of likely features down to points.
        [point_coordinates(:,2),point_coordinates(:,1)] = find(search_record) ; % Extract the locations of the identified features.
        %% END FEATURE IDENTIFICATION CODE BLOCK.
        
        if abort == 1
            clear point_coordinates
            clear abort
            break ;
        end
        
        % Records the total number of peaks found as a function of test box size.
        total_atoms = size(point_coordinates,1) ;
        k(1,((trialSize-1)/2)) = total_atoms ;
        
        k_diff  = gradient(k) ;
        
        % New automised ending approach. Note these are only active from the second
        % search and onwards, hence the first 'if' statement:
        
        if trialSize >= 15 % Value of 15 gives the minimum 5 unique points for quartic fitting.
            %Find quartic  fit:
            fittedGradient = polyfit( size_axis(1:(trialSize-5)/2) , log( -k_diff(1:((trialSize-5)/2)) ),4) ;
            [gradientLock,Index] = min( polyval( fittedGradient , size_axis( 1 : ((trialSize-5)/2) ) ) ) ;
            
            % Next some conditional statements to end the search.
            if ((trialSize-1)/2) > 1.25 * Index && log(-k_diff(1,((trialSize-5)/2))) > gradientLock*log(10)
                msgtitle = horzcat( 'Optimum feature spacing determined at ' , num2str(size_axis(1,Index)) , 'px. Total number of atoms is ' , num2str(k(Index)) , '.'  ) ;
                msgbox(msgtitle,'Ranger Results','help')
                clear point_coordinates
                subplot(1,2,2) % One final updating of the figure.
                semilogy(size_axis,k,'red')
                axis square
                xlim([0 big])
                ylim([10 numel(input_file)]/9)
                grid on
                xlabel ('Feature-Feature Spacing (px)','FontWeight','Bold')
                ylabel ('Normalised Amount','FontWeight','Bold')
                title('Number of Features Identified','FontWeight','Bold','FontSize',11)
                hold on
                plot(size_axis(1,Index),k(1,Index),'o','MarkerEdgeColor','blue')
                hold off
                pause(0.1)
                
                break % Search loop break.
            end
        elseif trialSize >= 13 % Value of 13 gives the minimum 4 unique points for cubic fitting.
            %Find cubic  fit:
            fittedGradient = polyfit( size_axis(1:(trialSize-5)/2) , log( -k_diff(1:((trialSize-5)/2)) ),3) ;
            [gradientLock,Index] = min( polyval( fittedGradient , size_axis( 1 : ((trialSize-5)/2) ) ) ) ;
            
            % Next some conditional statements to end the search.
            if ((trialSize-1)/2) > 1.25 * Index && log(-k_diff(1,((trialSize-5)/2))) > gradientLock*log(10)
                msgtitle = horzcat( 'Optimum feature spacing determined at ' , num2str(size_axis(1,Index)) , 'px. Total number of atoms is ' , num2str(k(Index)) , '.'  ) ;
                msgbox(msgtitle,'Ranger Results','help')
                clear point_coordinates
                subplot(1,2,2) % One final updating of the figure.
                semilogy(size_axis,k,'red')
                axis square
                xlim([0 big])
                ylim([10 numel(input_file)]/9)
                grid on
                xlabel ('Feature-Feature Spacing (px)','FontWeight','Bold')
                ylabel ('Normalised Amount','FontWeight','Bold')
                title('Number of Features Identified','FontWeight','Bold','FontSize',11)
                hold on
                plot(size_axis(1,Index),k(1,Index),'o','MarkerEdgeColor','blue')
                hold off
                pause(0.1)
                
                break % Search loop break.
            end
        elseif trialSize >= 11 % Value of 11 gives the minimum 3 unique points for quadratic fitting.
            %Find quadratic  fit:
            fittedGradient = polyfit( size_axis(1:(trialSize-5)/2) , log( -k_diff(1:((trialSize-5)/2)) ),2) ;
            [gradientLock,Index] = min( polyval( fittedGradient , size_axis( 1 : ((trialSize-5)/2) ) ) ) ;
            
            % Next some conditional statements to end the search.
            if ((trialSize-1)/2) > 1.25 * Index && log(-k_diff(1,((trialSize-5)/2))) > gradientLock*log(10)
                msgtitle = horzcat( 'Optimum feature spacing determined at ' , num2str(size_axis(1,Index)) , 'px. Total number of atoms is ' , num2str(k(Index)) , '.'  ) ;
                msgbox(msgtitle,'Ranger Results','help')
                clear point_coordinates
                
                subplot(1,2,2) % One final updating of the figure.
                semilogy(size_axis,k,'red')
                axis square
                xlim([0 big])
                ylim([10 numel(input_file)]/9)
                grid on
                xlabel ('Feature-Feature Spacing (px)','FontWeight','Bold')
                ylabel ('Number of Features','FontWeight','Bold')
                title('Number of Features Identified','FontWeight','Bold','FontSize',11)
                hold on
                plot(size_axis(1,Index),k(1,Index),'o','MarkerEdgeColor','blue')
                hold off
                pause(0.1)
                
                break % Search loop break.
            end
        else % This 'else' statement allows a gradient based match to be calcualted when insufficient uniquew points exist for quadratic fitting method above, but does not have the ability to teminate the loop.
            [~,Index] = max(k_diff) ;
        end
        
        %-----------------------------------------
        
        if show_progress == 1
            
            subplot(1,2,1)
            imagesc(input_file)
            colormap gray
            axis equal
            axis tight
            title('Coarse Positions','FontWeight','Bold','FontSize',11)
            set(gca,'XTick',[],'YTick',[])
            hold on
            plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
                'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
            hold off
            
            subplot(1,2,2)
            semilogy(size_axis,k,'red')
            axis square
            xlim([0 big])
            ylim([10 numel(input_file)]/9)
            grid on
            xlabel ('Feature-Feature Spacing (px)','FontWeight','Bold')
            ylabel ('Normalised Amount','FontWeight','Bold')
            title('Number of Features Identified','FontWeight','Bold','FontSize',11)
            hold on
            plot(size_axis(1,Index),k(1,Index),'o','MarkerEdgeColor','blue')
            hold off
            
            pause(0.1)
            
        end
        
        clear point_coordinates
        
        
    end
    
    delete(h2)
    % Once the optimum box size is found then this is applied once more.
    
    best_size = size_axis(1,Index) ;
    clear size_axis
    clear I
    
    if best_size == big
        msgtitle = horzcat( 'Upper test-box size limit reached. Size is ' , num2str(best_size) , 'px. Total number of atoms is ' , num2str(total_atoms) , '. Consider increasing upper limit to verify this result.'  ) ;
        msgbox(msgtitle,'Ranger: Warning','warn')
    end
    
    %         close (progFig)
    
end

% Variable clean up:
clear i
clear j
clear k
clear remaining_points
clear intensity_diff
clear k_diff
clear mean_intensity
clear fittedGradient
clear gradientLock
clear Index
clear abort
clear progFig

trialSize = best_size ;

%% BEGIN CODE FOR FEATURE IDENTIFICATION - SHOULD REMAIN IDENTICAL TO COPY ABOVE.
inputOffset = single(input_file - min(input_file(:))) ;         % Removes 'DC offset' from image to simplify Gaussian fitting.
vert_offset = zeros(m,n) ;                              % Create blank arrays.
horz_offset = zeros(m,n) ;
peak        = zeros(m,n) ;
spread      = zeros(m,n) ;
test_box_padding = ( trialSize - 1 ) / 2 ;              % Half of the trial size, equivalent to the border that will not be inspected.
baseAxis = ( -test_box_padding : test_box_padding )  ;  % Coordinate set for X and Y fitting.
baseAxis = baseAxis(:) ;
A = [baseAxis.^2 , baseAxis , ones(size(baseAxis))] ;
% Followed by the restoration progress bar:
h2=waitbar(0,'Identifying Image Peaks...','Name',version_string);
hw2=findobj(h2,'Type','Patch');
set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % Changes the color to red.
for i = test_box_padding + 1 : m - ( test_box_padding + 1 )
    
    currentStrip = inputOffset( i - test_box_padding : i + test_box_padding , : ) ;
    
    for j = test_box_padding + 1 : n - ( test_box_padding + 1 )
        I = currentStrip( : ,  j - test_box_padding : j + test_box_padding ) ;
        
        horizontalProfile = sum(I,1) ; % Generate two 1D profiles for Gaussian fitting.
        b = log(horizontalProfile(:)); % Natural log of the 1D data
        solutionH = A \ b ;                 % Least-squares solution
        horz_offset(i,j) = -solutionH(2)/solutionH(1)/2 ; % The horizontal offset refinement.
        
        verticalProfile = sum(I,2)' ; % Generate two 1D profiles for Gaussian fitting.
        b = log(verticalProfile(:)); % Natural log of the 1D data
        solutionV = A \ b ;                 % Least-squares solution
        vert_offset(i,j) = -solutionV(2)/solutionV(1)/2 ; % The vertical offset refinement.
        
        peak(i,j) = max(horizontalProfile(:)+verticalProfile(:))/(2*trialSize) ; % Calculates the height of the fitted Gaussian.
        spread(i,j) = real(sqrt( -1/2/(0.5*solutionH(1)+0.5*solutionV(1)) )) ; % Calculates the width of the fitted Gaussian.
    end
    
    percentageRefined = (i-test_box_padding) ./ (m - 2*test_box_padding) ;
    waitbar(percentageRefined,h2) ;
end

clear A
clear b
clear solutionH
clear solutionV
delete (h2)
% Feature identification section:
normalisedPeak = (peak ./ ( max(max(inputOffset)) - min(min(inputOffset))) ) ; % Make use of peak knowledge:
normalisedPeak(normalisedPeak<0) = 0 ; % Forbid negative (concave) Gaussians.
spread = spread ./ trialSize ;         % Normalise fitted Gaussian widths.
offsetRadius = realsqrt( (horz_offset).^2 + (vert_offset).^2 ) ; % Calculate offset radii.
offsetRadius = offsetRadius ./ trialSize ;
offsetRadius(offsetRadius==0) = 0.001 ; % Remove zeros values to prevent division error later.
% Create search metric and screen impossible peaks:
search_record = ( normalisedPeak ./ offsetRadius ) ;
search_record(search_record > 1) = 1 ;
search_record(search_record < 0) = 0 ;
search_record(spread < 0.05) = 0 ;      % Invalidates negative Gaussian widths.
search_record(spread > 1) = 0 ;         % Invalidates Gaussian widths greater than a feature spacing.
search_record(offsetRadius > 1) = 0 ;   % Invalidates Gaussian widths greater than a feature spacing.
search_record = medfilt2(search_record, [round(trialSize/3) round(trialSize/3)]) ; % Median filter to strip impossibly local false-positive features.
sensitivityThreshold = 0.34 ;           % This is an Admin tunable parameter that is defined here within the core file.
search_record(search_record < sensitivityThreshold ) = 0 ;  % Collapse improbable features to zero likelyhood.
search_record(search_record >= sensitivityThreshold ) = 1 ; % Round likelyhood of genuine features to unity.
% Logical matrix is stored at this point for horizontal calculation acceleration...
informationRegion = search_record ;
% ...after which search_record operations can continue.
search_record = bwmorph(search_record,'shrink',Inf) ;       % Errod regions of likely features down to points.
[point_coordinates(:,2),point_coordinates(:,1)] = find(search_record) ; % Extract the locations of the identified features.
%% END FEATURE IDENTIFICATION CODE BLOCK.

% Variable cleanup:
clear inputOffset
clear peak
clear spread
clear normalisedPeak
clear offsetRadius
clear currentStrip
clear sensitivityThreshold

%% Manual Exception Section.
% The following sub-section allows for extra peaks to be manually added.

% To begin with we will need a figure to view the image:

manSelec = figure('Position',[58 148 fullscreen(3)-115 fullscreen(4)-336],'Name',version_string,'NumberTitle','off' ) ;

imagesc(inputUnfiltered)
axis image
colormap gray
title ('Computer found peaks - click any missing peaks.','FontWeight','Bold','FontSize',11)
hold on
plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
    'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
hold off
set(gca,'XTick',[],'YTick',[])

still_editing = 1 ; % This will determine whether the user is still adding or removing points.

while still_editing == 1
    
    [a,b] = getpts(gcf) ; % Gets the extra peaks needed from the cursor clicks.
    extras = round(cat(2, a, b)) ; % Merge the x and y coordinates of the extra peaks...
    clear a
    clear b
    
    for extra = 1:size(extras,1) % Update the search-record for forward comaptability.
        search_record(extras(extra,2),extras(extra,1)) = 1 ;
    end
    
    point_coordinates = cat(1, point_coordinates, extras) ; % ... and add these onto the list of coordinates.
    
    % Update the figure:
    pause(0.25)
    
    imagesc(inputUnfiltered)
    axis image
    colormap gray
    title ('Computer found peaks - click any excess peaks.','FontWeight','Bold','FontSize',11)
    hold on
    plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
        'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
    hold off
    set(gca,'XTick',[],'YTick',[])
    
    % The following sub-section allows peaks to be manually removed.
    
    [a,b] = getpts(gcf) ; % Gets the extra peaks needed from the cursor clicks.
    surplus = round(cat(2, a, b)) ; % Merge the x and y coordinates of the extra peaks...
    
    % Now we have a list of the peaks that we'd like to remove. They may
    % however not exactly match the coordinates of pre-exisiting peaks.
    % Instead we will remove the nearest pre-exisiting peak to those
    % identified.
    for scrub = 1: size(surplus,1)
        % First calculate the distance to all the pre-existing peaks:
        distance = realsqrt((point_coordinates(:,1)-surplus(scrub,1)).^2 + (point_coordinates(:,2)-surplus(scrub,2)).^2) ;
        % Next find the index that corresponds to that closest point:
        [~,I] = min(distance) ;
        % and strip that coordinate out:
        search_record(round(point_coordinates(I,2)),round(point_coordinates(I,1))) = 0 ; % Update the search-record for forward comaptability.
        point_coordinates(I,:) = [] ;
        clear distance
    end
    
    % Update the figure:
    pause(0.25)
    
    imagesc(inputUnfiltered)
    axis image
    colormap gray
    title ('Computer found peaks - click any missing peaks.','FontWeight','Bold','FontSize',11)
    hold on
    plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
        'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
    hold off
    set(gca,'XTick',[],'YTick',[])
    
    if size(extras,1) == 0 && size(surplus,1) == 0
        still_editing = 0 ;
    end
    
end

title ('Computer found peaks - end result.','FontWeight','Bold','FontSize',11)

pause(1)
close(manSelec)
% Clean up:
clear manSelec
clear a
clear b
clear extra
clear extras
clear scrub
clear surplus
clear still_editing

% END Manual Exception Section
%%


% Update the search record with the new points:
search_record = 0 * search_record ; % Erases the current search_record without changing its dimensions.
for counter = 1:size(point_coordinates,1)
    search_record(round(point_coordinates(counter,2)),round(point_coordinates(counter,1))) = inputUnfiltered(round(point_coordinates(counter,2)),round(point_coordinates(counter,1))) ;
end
clear counter
search_record = 1 * search_record ; % Converts the logical array search record to ones and zeros.



clear trialSize
clear i
clear j
clear m
clear n
clear horz_offset
clear vert_offset
clear test_pix
clear test_box_padding
clear test_box_size
clear percentageRefined

total_atoms = size(point_coordinates,1) ;

% BEGIN Position Refinement Section

if refinePositions == 1
    % This section generates the coordinates vector for the atom positions.
    [row,col] = find(search_record) ;
    total_atoms = size(row,1) ; % Updates number of total peaks in the case there are duplicates.
    coordinates(:,1) = col ;
    coordinates(:,2) = row ;
    % Check that no points fall too close to the image edges:
    col(col < (best_size+1)/2) = (best_size+1)/2 ;
    col(col > size(inputUnfiltered,2)-(best_size+1)/2) = size(inputUnfiltered,2)-(best_size+1)/2 ;
    row(row < (best_size+1)/2) = (best_size+1)/2 ;
    row(row > size(inputUnfiltered,1)-(best_size+1)/2) = size(inputUnfiltered,1)-(best_size+1)/2 ;
    % Create blank vectors:
    horz_offset = zeros(1,total_atoms) ;
    vert_offset = zeros(1,total_atoms) ;
    col_refined = zeros(1,total_atoms) ;
    row_refined = zeros(1,total_atoms) ;
    
    % Followed by the restoration progress bar:
    h2=waitbar(0,{'Now Performing Position Refinement.'},'Name',version_string);
    hw2=findobj(h2,'Type','Patch');
    set(hw2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0]) % changes the color to red.
    
    baseAxis = ( -(best_size-1)/2 : (best_size-1)/2 )  ; % Coordinate set for X and Y fitting.
    baseAxis = baseAxis(:) ;
    A = [baseAxis.^2 , baseAxis , ones(size(baseAxis))];
    
    for feature = 1 : size(row,1)
        I = inputUnfiltered( (row(feature) - (best_size-1)/2) : (row(feature) + (best_size-1)/2) , ...
            (col(feature) - (best_size-1)/2) : (col(feature) + (best_size-1)/2) ) ;
        
        horizontalProfile = sum(I,1) ; % Generate two 1D profiles for Gaussian fitting.
        b = log(horizontalProfile(:)); % Natural log of the 1D data
        solution = A \ b ;                 % Least-squares solution for x
        if abs(-solution(2)/solution(1)/2) < (best_size/3) % These lines check whether the suggested offsets are physically real:
            horz_offset(feature) = -real(solution(2)/solution(1)/2) ;
        else
            horz_offset(feature) = 0 ;
        end
        
        verticalProfile = sum(I,2)' ; % Generate two 1D profiles for Gaussian fitting.
        b = log(verticalProfile(:)); % Natural log of the 1D data
        solution = A \ b ;                 % Least-squares solution for x
        if abs(-solution(2)/solution(1)/2) < (best_size/3) % These lines check whether the suggested offsets are physically real:
            vert_offset(feature) = -real(solution(2)/solution(1)/2) ;
        else
            vert_offset(feature) = 0 ;
        end
        
        col_refined(feature) = col(feature) + horz_offset(feature) ; % Update the coordinate.
        row_refined(feature) = row(feature) + vert_offset(feature) ; % Update the coordinate.
        
        
        if (feature/10) == round(feature./10)
            percentageRefined = feature ./ size(coordinates,1) ;
            waitbar_progress_text1 = horzcat ( 'Now refining position ' , num2str(feature) , ' of ' , num2str(size(coordinates,1)) , '.' ) ;
            waitbar_progress_text2 = horzcat ( 'Refinement progress = ' , num2str(100*percentageRefined) , '%.' ) ;
            waitbar(percentageRefined,h2,{waitbar_progress_text1;waitbar_progress_text2},'Name',version_string) ;
        end
    end
    
    % Variable clean up:
    clear row
    clear col
    clear coordinates
    clear I
    clear m
    clear n
    clear feature
    clear A
    clear b
    clear baseAxis
    clear horizontalProfile
    clear verticalProfile
    clear waitbar_progress_text1
    clear waitbar_progress_text2
    clear percentageRefined
    clear solution
    
    figure
    plot(horz_offset, vert_offset,'o', ...
        'MarkerSize', 6, 'MarkerEdgeColor','k', 'MarkerFaceColor','b')
    grid on
    axis square
    xlabel('Horizontal Offset','FontWeight','Bold')
    ylabel('Vertical Offset','FontWeight','Bold')
    
    delete (h2)
    clear horz_offset
    clear vert_offset
    
    % For compatability with other modules an integer pixel coordinate
    % search_record should be returned:
    search_record = 0 * search_record ; % Erases the current search_record without changing its dimensions.
    for counter = 1:size(row_refined,2)
        search_record(round(row_refined(counter)),round(col_refined(counter))) = inputUnfiltered(round(row_refined(counter)),round(col_refined(counter))) ;
    end
    clear counter
end

search_record = single(search_record) ; % Converts data type to single precision.

% END Position Refinement Section

if display_results == 1
    
    figure('Position',[58 148 fullscreen(3)-115 fullscreen(4)-336],'Name',version_string,'NumberTitle','off' ) ;
    
    imagesc(inputUnfiltered)
    colormap gray
    axis equal
    axis tight
    title('Coarse Positions','FontWeight','Bold','FontSize',11)
    set(gca,'XTick',[],'YTick',[])
    hold on
    plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
        'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
    hold off
    
    % Create new figure and overlay points on image:
    
    if refinePositions == 1
        clear point_coordinates
        point_coordinates(:,1) = col_refined ;
        point_coordinates(:,2) = row_refined ;
        clear col_refined
        clear row_refined
        
        figure('Position',[58 148 fullscreen(3)-115 fullscreen(4)-336],'Name',version_string,'NumberTitle','off' ) ;
        imagesc(inputUnfiltered)
        axis image
        colormap gray
        title ('Refined Positions','FontWeight','Bold','FontSize',11)
        hold on
        plot(point_coordinates(:,1),point_coordinates(:,2),'o', ...
            'MarkerSize', 2, 'MarkerEdgeColor','k', 'MarkerFaceColor',[.49 1 .63])
        hold off
        set(gca,'XTick',[],'YTick',[])
    end
    
    if best_size < big
        msgtitle = horzcat( 'Optimum box size found. Size is ' , num2str(best_size) , 'px. Total number of atoms is ' , num2str(total_atoms) , '.'  ) ;
        rang_box = msgbox(msgtitle,version_string,'help') ;
        pause(1)
        close(rang_box)
        clear rang_box
    end
    
    clear ranger_progress
end

clear msgtitle
clear test_area
clear big
clear h2
clear hw2

% Before closing this scrip any image filtering should be overwritten
% and the filtered image discarded:
input_file = inputUnfiltered ;
clear inputUnfiltered
