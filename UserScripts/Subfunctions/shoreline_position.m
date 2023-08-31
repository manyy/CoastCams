function [shoreline_positions] = shoreline_position (I, threshold, res, ShoreMethod, plotoption_shore)
% Inputs:
% I = Timestack Image
% ShoreMethod:
%   1 = grayscale method specifically designed for the rocky platform
%   2 = Red- minus Blue-channel method
%   3 = Colour divergence method

%I = imread('G:\CAMCOAST-master\CAMS\CAMS_DATA\02 - DATA\GRANDPOPO\GPP_NIVEAU 1\S_3_202108041000.jpg');
%ShoreMethod = 3;

if ShoreMethod == 1; % Grayscale method specifically for rocky platform Socoa
    gray_av = [];
    gray = rgb2gray(I);
    gray_adj = imadjust(gray, stretchlim(gray), [0 1], 1);
    gray_av = [gray_av; mean(vertcat(gray_adj),1,'omitnan')];
    J = imrotate (gray_av, 90);
    
    threshold = 30;
    
    [r c q]=size(J);
    
    xti = [1 size(J,2) size(J,2) 1]; 
    yti = [0, 0, ceil(size(J,1) ),ceil(size(J,1))]; 
        
    P = improfile(J, xti, yti);
    F = flip(P);
    Ps = movmean(F, 20);
    m = 1:length(Ps);
               
    [mins, min_locs] = findpeaks(-Ps, 'MinPeakProminence', 8);
    [maxs, max_locs] = findpeaks(Ps, 'MinPeakProminence', 8);
    
    peaks = [0 0; maxs max_locs];
    troughs = [mins min_locs];
    
    sp = size (peaks);
    st = size (troughs);
    s = max(sp(1), st(1));
    dif = [[peaks;zeros(abs([s 0] - sp))],[troughs;zeros(abs([s, 0]-st))]];
    difs = [dif (dif(:,1)+dif(:,3))]; 
             
    S10 = imrotate(J, 180);     
       
    idx = find(difs(:,5)>threshold,1);
    if ~isempty(idx)
        shoreline_positions = difs(idx,4);
    else shoreline_positions = 0;
    end  
    
    if plotoption_shore == 1;
        f1=figure;
        imagesc(I)
        yyaxis right
        hold on
        plot (m, Ps, 'w')
        plot (min_locs, Ps(min_locs), 'rv', 'MarkerFaceColor', 'r')
     	plot(max_locs, Ps(max_locs), 'rs', 'MarkerFaceColor', 'b')
        %plot(shore, Ps(shore), 'rv', 'MarkerFaceColor', 'g');
        set(gca, 'fontsize',14)
    end               
end

if ShoreMethod == 2; %Color Channel Divergence
    ccd_av = [];
    ccd_av = [ccd_av; mean(vertcat(I),1,'omitnan')];
    J_ccd = imrotate (ccd_av, 90); 
    
    xtcc = [1 size(J_ccd,2) size(J_ccd,2) 1]; 
    ytcc = [0, 0, ceil(size(J_ccd,1) ),ceil(size(J_ccd,1))]; 
        
    [pdf_values,pdf_locs] = ksdensity(J_ccd(:,:,1)-J_ccd(:,:,3)); %find smooth pdf of CCD (Red minus Blue) channel
    xlabel_type = 'Red minus blue';
    %thresh_weightings = [1/3 2/3]; %This weights the authomatic thresholding towards the sand peak (works well in SE Australia)
    thresh_weightings = [1/3 2/3];
    %[peak_values,peak_locations]=findpeaks(pdf_values,pdf_locs); %Find peaks
    thresh_otsu = multithresh(J_ccd(:,:,1)-J_ccd(:,:,3)); %Threshold using Otsu's method
    I1 = find(pdf_locs<thresh_otsu);
    [~,J1] = max(pdf_values(I1));
    I2 = find(pdf_locs>thresh_otsu);
    [~,J2] = max(pdf_values(I2));
    thresh = thresh_weightings(1)*pdf_locs(I1(J1)) + thresh_weightings(2)*pdf_locs(I2(J2)); %Skew average towards the positive (i.e. sand pixels)
    disp(['Sand/water threshold determined to be ' num2str(thresh, '%0.0f')])
    
    if plotoption_shore==1
        f2 = figure;
        plot(pdf_locs,pdf_values)
        hold on
        plot(pdf_locs([I1(J1) I2(J2)]),pdf_values([I1(J1) I2(J2)]),'ro')
        YL = ylim;
        plot([thresh thresh], YL,'r:','linewidth',2)
        %plot([thresh_otsu thresh_otsu], YL,'g:','linewidth',2)
        xlabel(xlabel_type,'fontsize',10)
        ylabel('Counts','fontsize',10)
    end
    
    %Extract contour
    RmB = double(I(:,:,1))- double(I(:,:,3));
    [C, h] = contours(RmB, [thresh thresh]);
    II = find(C(1,:)==thresh);
    
if ShoreMethod == 3
    % Convert the timestack image to a grayscale image
gray_I = rgb2gray(I);

% Calculate the red minus blue channel image
red_I = double(I(:,:,1));
blue_I = double(I(:,:,3));
RMB_I = red_I - blue_I;

% Threshold the RMB image to create a binary image
level = graythresh(RMB_I);
bw_RMB = imbinarize(RMB_I, level);

% Remove noise and small objects from the binary image
se = strel('disk', 5);
bw_RMB = imopen(bw_RMB, se);

% Identify the shoreline in each column of the binary image
shoreline = zeros(1, size(bw_RMB, 2));
for col = 1:size(bw_RMB, 2)
    row = find(bw_RMB(:,col), 1, 'first');
    if isempty(row)
        shoreline(col) = NaN;
    else
        shoreline(col) = row;
    end
end

% Convert the shoreline positions from pixel coordinates to physical coordinates
% (assuming a scale of 1 pixel = 1 meter)
shoreline_x = (1:size(bw_RMB, 2));
shoreline_y = shoreline;

% Plot the shoreline position on the timestack image
figure;
imshow(I);
hold on;
plot(shoreline_x, shoreline_y, 'r', 'LineWidth', 2);
end
    
  
end
    
    
    
    
    
    