% function to calculate correlation between two 3D array
% 
function [G2, G2_blur] = fn_calc_G2(data1, data2, norm1)

% reshape 3D matrix to 2D with: raws are pixel region, columns - time axis
if size(data1) > 2 
    data1 = reshape(data1, size(data1,1)*size(data1,2), size(data1,3));
end
if size(data2) > 2 
    data2 = reshape(data2, size(data2,1)*size(data2,2), size(data2,3));
end

% check in case there is nans and remove it; 
% todo: need to be corrected for correlation of different datasets: remove all nans_data1 * nans_data2
indnan = [];
for nn = 1:1:size(data1)
	if isnan(data1(nn,1))
	indnan = [indnan nn];
	end
end

data1(indnan,:) = [];
data2(indnan,:) = [];

% if normalization is set to 1
if norm1 == 1
    data1 = data1./squeeze(max(data1,[],1));
    data2 = data2./squeeze(max(data2,[],1));
end


% for parallelisation use another updated code
n =size(data1,2);
g2 = zeros(n,n);
g22 = [];

%parfor ii = 1:n
for ii = 1:n
    jj = ii-1;
    I1 = [];
    I2 = [];
    I1I2 = [];
    meanI1I2 = [];
    g2_help = [];

    I1 = data1(:,1:end-jj);
    I2 = circshift(data2, -jj, 2);
    I2 = I2(:,1:end-jj);
    
% I guess you want to add something here like: if jj == 0:


%     I1I2 = (I1).*(I2);
%     meanI1I2 = mean(I1I2,1);
%     %sqrtI1I2 = sqrt(mean((I1).^2,1).*mean((I2).^2,1));
%     g2_help = (meanI1I2);%./sqrtI1I2;
%     g2(ii,1+jj:end) = g2_help'; 
%     g22(ii) = mean(g2_help);

%type 2 
    meanI1rep = [];
    meanI2rep = [];    
    meanI1 = mean(I1,1);
    meanI2 = mean(I2,1);
    meanI1rep = repmat(meanI1, size(I1,1), 1);
    meanI2rep = repmat(meanI2, size(I2,1), 1);
 %%%%%%%%% std / sigma
    I1I2 = (I1-meanI1rep).*(I2-meanI2rep);
    sumI1I2 = mean(I1I2,1);
    sqrtI1I2 = sqrt(mean((I1-meanI1rep).^2,1).*mean((I2-meanI2rep).^2,1));
    g2_help = (sumI1I2)./sqrtI1I2;
    g2(ii,1+jj:end) = g2_help'; 
    g22(ii) = mean(g2_help);

    
%Eric's Dufresne way
%     I1I2 = ((I1-meanI1rep)./meanI1rep).*((I2-meanI2rep)./meanI2rep);
%     sumI1I2 = mean(I1I2,1);
%     %%%%%%%%%% sqrtI1I2 = sqrt(mean((I1-meanI1rep).^2,1).*mean((I2-meanI2rep).^2,1));
%     g2_help = sumI1I2;
%     g2(ii,1+jj:end) = g2_help'; 
%     g22(ii) = mean(g2_help);

end

% fill two-time matrix
G2 = zeros(size(g2));
for ii = 1:1:size(G2,1)
    for jj = 1:1:size(G2,2)
        dif = abs(ii-jj);    
        G2(ii,jj) = g2(dif+1,max(ii,jj));
        
    end
end

% change main diagonal
G2_blur = G2;
G2_blur(1,1) = mean([G2(2,1) G2(1,2)]);
G2_blur(size(G2,1),size(G2,2)) = mean([G2(size(G2,1)-1,size(G2,2)) G2(size(G2,1),size(G2,2)-1)]);
for ii = 2:1:size(G2,1)-1
    jj = ii;
    G2_blur(ii,jj) = mean([G2(ii+1,jj) G2(ii-1,jj) G2(ii,jj+1) G2(ii,jj-1)]);
end


% calculate g2 from G2_blur
differ_vector = 0:1:size(G2,1)-1;
row_help_map = repmat(1:1:size(G2,1),size(G2,1),1);
coloumn_help_map = transpose( repmat(1:1:size(G2,1),size(G2,1),1));
diag_map = abs(row_help_map-coloumn_help_map);

for ii = 1:1:size(G2,1)%length(diff_vector)
    diag_value = ii-1;
    g2_from_G2blur(ii) = mean(G2_blur(diag_map==diag_value));
    
end


% to make plots
% timepreset = 0.25;
% 
% timediag = diag_map*timepreset;
% y = 1:1:size(data1,3);
%     y = y-1;
%     y = y*timepreset;
% % logtimescale = 
% 
% disp('done')

% %%
% figure(20)
% clf
% set(gcf, 'color', 'white', 'Position', [526   579   515   448])
% subplot(2,2,1)
% imagesc(mean(data1,3));
% xlabel('pixel')
% ylabel('pixel')
% hold all
% axis image
% set(gca, 'fontsize', 14)
% colorbar
% title('ROI')
% %plot([ROI_3 ROI_3 ROI_4 ROI_4 ROI_3], [ROI_1 ROI_2 ROI_2 ROI_1 ROI_1], 'r')
% %imagesc(data(80:180,108:149,1))
% subplot(2,2,2)
% imagesc(G2)
% colorbar
% axis image
% axis xy
% xlabel('frame')
% ylabel('frame')
% set(gca, 'fontsize', 14)
% 
% subplot(2,2,3)
% imagesc(G2_blur)
% colorbar
% axis image
% axis xy
% xlabel('frame')
% ylabel('frame')
% set(gca, 'fontsize', 14)
% 
% 
% subplot(2,2,4)
% plot(y, g2_from_G2blur, '-', 'linewidth', 1);
% set(gca, 'Xscale', 'log')
% xlabel('frame')
% ylabel('g2')
% %ylim([0.03, 0.05])
% set(gca, 'fontsize', 14)
% 
% 
% 
% set(gcf, 'Position', [ 567          96        1115         871])
% 
% 
% %% take delta perpendicular to main diagonal and convert it to linear
% steps plot

% figure()
% deltaNframe = 100;
% diag_map_help = zeros(size(diag_map));
% diag_map_help(diag_map<deltaNframe) =1;
% G2_blur_nonzeros = G2_blur;
% G2_blur_nonzeros(G2_blur_nonzeros<0) = 0;
% imagesc(G2_blur_nonzeros.*diag_map_help)
% colorbar
% axis image
% axis xy
% xlabel('frame')
% % histogramming
% 
% G2_hist_region = G2_blur(diag_map_help==1);
% 
% figure()
% hist(G2_hist_region)
% %
% linear_two_time_represantation = zeros(deltaNframe+1,size(G2,1));
% rectsize1 = size(size(G2,1));
% rectsize2 = deltaNframe+1;
% for ii = 1:1:deltaNframe+1;
%     find(diag_map ==ii-1);
%     [i_ind,j_ind] =  ind2sub(size(G2),find(diag_map ==ii-1));
%     diffxy = i_ind-j_ind;
%     pos_helpi = i_ind(diffxy>=0);
%     pos_helpj = j_ind(diffxy>=0);
%     index_for_rect = sub2ind(size(G2), pos_helpi,pos_helpj);
%     
%     linear_two_time_represantation(ii,ii:(ii-1)+length(index_for_rect)) = G2_blur(index_for_rect);
%     
% end
% 
% figure
% imagesc(linear_two_time_represantation)
% axis xy
