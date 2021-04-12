clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mapping Sample Space to Detector Conjugated Space
% Version 1.0
% May 2019
% Written By: David Yang
% University of Oxford, Dept. of Engineering Science
% 
% PURPOSE: To map an object in sample space (SS) to detector conjugated space (DCS) for the same object.
% 
% FILE INPUT: A file containing a 3D sample space array
% 
% FILE OUTPUT: A file containing the calculated DCS in a -BIN.mat file
% 
% USER-DEFINED SECTIONS:
% 1. Navigating to the script folder
% - changes the MATLAB active folder to a user-specified directory
% 
% 2. Original reflection details
% - collects details about the reconstruction DCS shape to be translated
% 
% 3. Beamline selection
% - user specifies the beamline plugin to create specific rotation matrices
% 
% 4. Plot options
% - user can plot the calculated shape and configure parameters
% 
% 5. Save calculated shape
% - option to save the calculated DCS shape or twin in a -BIN.mat file
% 
% 6. Test
% - to plot the calculated shape with a previous reconstruction for comparison
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
fprintf('       Mapping Sample Space to Detector Conjugated Space\n');
fprintf('                         Version 1.0\n');
fprintf('                           May 2019\n');
fprintf('                    Written By: David Yang\n');
fprintf('      University of Oxford, Dept. of Engineering Science\n');
fprintf('<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>---<>\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Original reflection details
fprintf('\n...collecting original reflection details...');

% directory to the reconstruction folder
O.dir = 'Examples/Cylinder_(-120)_-36.0786_gamma_40.2954_delta_0.00274_dtheta';

% name of the reconstruction folder
O.file_name = 'Cylinder_(-120)_-36.0786_gamma_40.2954_delta_0.00274_dtheta-SAM';

% importing the SS shape
try
    O.SS_shape_REC = cell2mat(struct2cell(load([O.dir, '/', O.file_name, '.mat'])));
catch
    O.SS_shape_REC = cell2mat(struct2cell(load([O.file_name, '.mat'])));
end

% the old hkl reflection
O.hkl = [-1; 2; 0];

% wavelength in m for original reflection
O.lambda = (12.398/10.0)/10*10^-9;

% get size of the original matrix
[O.N1, O.N2, O.N3] = size(O.SS_shape_REC);
O.p_sam = 1; % obtained from reconstruction, set to 1 if unknown

% beamline sample motor angles in degrees (set to the motors' default angles for lab space)
% O.theta_bl = 0; % for lab space
% O.chi_bl = 90; % for lab space
% O.phi_bl = 0; % for lab space
O.theta_bl = -25.9718;
O.chi_bl = 153.4349;
O.phi_bl = -47.0851;

% beamline detector motor angles in degrees
O.delta_bl = 40.2954;
O.gamma_bl = 36.0786;

% choose rocking angle and increment in degrees
O.rocking_angle = 'dtheta'; % 'dphi' rotate about x-axis, 'dtheta' rotate about y-axis for 34-ID-C
O.rocking_increment = 0.00274; % rocking angle step size

% detector parameters in m
O.D = 1.7745; % detector distance in m, put absolute value if known
O.d = 5.5e-5; % detector pixel size in m--this is effective, becomes larger than actual if binning is used
O.N = 256; % number of pixels along one dimension of square detector


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Beamline selection
% adding path to APS 34-ID-C angles folder
addpath(genpath('APS 34-ID-C angles'));

% beamline-specific plugin for O
[O.R_dqp_12, O.R_dqp_3, O.R_xyz, O.S_0lab_dir] = plugin_APS_34IDC(O.theta_bl, O.chi_bl, O.phi_bl, O.delta_bl, O.gamma_bl, O.rocking_increment, O.rocking_angle);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Plot options
% plot the calculated DCS shape
plot_DCS_shape = 1; % 1 to plot calculated DCS shape

% toggle between viewing the thresholded amplitude or isosurface
amplitudes = 1; % 1 to plot amplitudes, 0 to plot isosurfaces
amplitude_threshold = 0.3; % 0.3 is a good starting point

% select viewpoint
viewpoint = [-180, -90]; % viewpoint = [az, el], x-y plane
% viewpoint = [-180, 0]; % viewpoint = [az, el], x-z plane
% viewpoint = [-90, 0]; % viewpoint = [az, el], y-z plane


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Save calculated shape
% choose to save the reflection
save_reflection = 1; % 1 to save, 0 to not save

% choose to take the twin of the new calculated DCS shape
twin = 1; % 1 to take the conjugate reflection

% choose to binarize the calculated DCS shape based on a threshold value
binarize = 1; % 1 to binarize according to binarize_threshold
binarize_threshold = 0.3;

% save directory
save_dir = 'Simulated Data'; 

% save name
save_name = 'Cylinder'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Testing
 % 1 to test against a reconstructed shape, 0 to not test
test = 1;

% directory to the reconstruction folder
O.dir = 'Examples';

% name of the reconstruction folder
O.file_name = 'Rec-Cylinder_(-120)_dtheta-00274-ERlrHIOlr2000-NM-SW';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% DON'T EDIT BELOW THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Making grids assuming 34-ID-C coordinate frame (DO NOT TOUCH)
fprintf('\n...making pixel grids...');
% make pixel coordinate grids for the 3D volume of the original DCS shape
[O.N1grid, O.N2grid, O.N3grid] = meshgrid(-(O.N1-1)/2:(O.N1-1)/2, -(O.N2-1)/2:(O.N2-1)/2, -(O.N3-1)/2:(O.N3-1)/2);

% calculating original voxel size if unknown
if O.p_sam == 1
    O.p_sam = O.lambda*O.D/(O.N*O.d);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make original (O) coordinates (DO NOT TOUCH)
fprintf('\n...making original detector conjugated space coordinates...');
% making S_0lab, S_lab and Q_lab vectors for beamline
O.S_0lab = 2*pi/O.lambda*O.S_0lab_dir;
O.S_lab = O.R_dqp_12*O.S_0lab;
O.Q_lab = O.S_lab - O.S_0lab;

% making x_sam , y_sam , z_sam sample space vectors
O.x_sam = O.R_xyz*[1; 0; 0];
O.y_sam = O.R_xyz*[0; 1; 0];
O.z_sam = O.R_xyz*[0; 0; 1];

% making q_1', q_2', q_3' detector reciprocal space vectors
O.q_1p = 2*pi/O.lambda*O.d/O.D*O.R_dqp_12*[1; 0; 0];
O.q_2p = 2*pi/O.lambda*O.d/O.D*O.R_dqp_12*[0; 1; 0];
O.q_3p = O.R_dqp_3*O.Q_lab-O.Q_lab;

% make x', y', and z' detector conjugated space basis vectors, adapted from Berenguer et al. PRB 88, 144101 (2013).
O.V_DRS = dot(cross(O.q_3p, O.q_2p), O.q_1p)*O.N1*O.p_sam*O.N2*O.p_sam*O.N3*O.p_sam;
O.xp = 2*pi*cross(O.N2*O.p_sam*O.q_2p, O.N3*O.p_sam*O.q_3p)./O.V_DRS;
O.yp = 2*pi*cross(O.N3*O.p_sam*O.q_3p, O.N1*O.p_sam*O.q_1p)./O.V_DRS;
O.zp = 2*pi*cross(O.N1*O.p_sam*O.q_1p, O.N2*O.p_sam*O.q_2p)./O.V_DRS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Map from original DCS to new DCS (DO NOT TOUCH)
fprintf('\n...interpolating sample space to detector conjugated coordinates...');
% making the T_DCS_to_SS matrix for O
O.T_DCS_to_SS = [O.x_sam, O.y_sam, O.z_sam]\[O.xp, O.yp, O.zp]; % equivalent to [a_1, a_2, a_3]*inv(b_1, b_2, b_3])

% map original coordinates to non-orthogonal coordinates
O.N1gridp = O.T_DCS_to_SS(1,1)*O.N1grid + O.T_DCS_to_SS(1,2)*O.N2grid + O.T_DCS_to_SS(1,3)*O.N3grid;
O.N2gridp = O.T_DCS_to_SS(2,1)*O.N1grid + O.T_DCS_to_SS(2,2)*O.N2grid + O.T_DCS_to_SS(2,3)*O.N3grid;
O.N3gridp = O.T_DCS_to_SS(3,1)*O.N1grid + O.T_DCS_to_SS(3,2)*O.N2grid + O.T_DCS_to_SS(3,3)*O.N3grid;

% interpolate sample space data to detector conjugated frame
O.DCS_shape_CALC = interp3(O.N1grid, O.N2grid, O.N3grid, O.SS_shape_REC, O.N1gridp, O.N2gridp, O.N3gridp, 'linear',  0); % make any values outside original data zero.
  
% taking the twin if necessary
if twin == 1
    F = ifftshift(fftn(fftshift(O.DCS_shape_CALC)));
    O.DCS_shape_CALC = fftshift(ifftn(ifftshift(conj(F))));
end

% centring about the centre of mass
O.DCS_shape_CALC_MASK = single(abs(O.DCS_shape_CALC) > binarize_threshold);
structure_element = strel('sphere', 3);
O.DCS_shape_CALC_MASK = imerode(imdilate(O.DCS_shape_CALC_MASK, structure_element),structure_element); % takes care of dislocation cores
O.DCS_shape_CALC_COM = ceil(centerOfMass(O.DCS_shape_CALC_MASK));
O.DCS_shape_CALC = circshift(O.DCS_shape_CALC, size(O.DCS_shape_CALC)/2-O.DCS_shape_CALC_COM);

% binarizing the calculated DCS shape
if binarize == 1
    fprintf(['\n...binarizing new reflection DCS shape with a threshold of ', num2str(binarize_threshold), '...']);
    O.DCS_shape_CALC = single(abs(O.DCS_shape_CALC) > binarize_threshold); % binarizing
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Saving new reflection shape (DO NOT TOUCH)
if save_reflection == 1
    fprintf('\n...saving new reflection shape...');
    mkdir(save_dir);
    array = O.DCS_shape_CALC;
    save([save_dir, '/', save_name,'_(', num2str(O.hkl(1)), num2str(O.hkl(2)), num2str(O.hkl(3)),')_', num2str(O.gamma_bl),'_gamma_',num2str(O.delta_bl),'_delta_',num2str(O.rocking_increment),'_',O.rocking_angle, '-BIN.mat'], 'array');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Testing calculated shape with reconstructed shape (DO NOT TOUCH)
if plot_DCS_shape == 1
    if test == 1
        fprintf('\n...plotting reconstructed detector conjugated space shape from reconstruction to compare...');
        
        % importing the -AMP.mat and -PH.mat files
        try
            O.DCS_shape_REC_AMP = cell2mat(struct2cell(load([O.dir, '/', O.file_name, '/', O.file_name,'-AMP.mat'])));
            O.DCS_shape_REC_PH = cell2mat(struct2cell(load([O.dir, '/', O.file_name, '/', O.file_name,'-PH.mat'])));
        catch
            fprintf('\n...the reconstruction folder and/or file is corrupt or cannot be found...');
            test = 0;
        end
        
        % generating the reconstructed DCS shape
        O.DCS_shape_REC = abs(O.DCS_shape_REC_AMP.*exp(1i.*O.DCS_shape_REC_PH));
        
        % choose to binarize the reconstructed DCS shape
        if binarize == 1
            O.DCS_shape_REC = single(abs(O.DCS_shape_REC) > binarize_threshold); % binarizing
        end
        
        % centring about the centre of mass
        fprintf('\n...centring DCS shapes...');
        O.DCS_shape_REC_MASK = single(abs(O.DCS_shape_REC) > amplitude_threshold);
        structure_element = strel('sphere', 3);
        O.DCS_shape_REC_MASK = imerode(imdilate(O.DCS_shape_REC_MASK, structure_element),structure_element); % takes care of dislocation cores
        O.DCS_shape_REC_COM = ceil(centerOfMass(O.DCS_shape_REC_MASK));
        O.DCS_shape_REC = circshift(O.DCS_shape_REC, size(O.DCS_shape_REC)/2-O.DCS_shape_REC_COM);
    end
    
    % creating a figure
    figure;
    hold on;
    
    if amplitudes == 1
        % plotting calculated shape amplitude
        O.plot = patch(isosurface(O.N1grid*O.p_sam, O.N2grid*O.p_sam, O.N3grid*O.p_sam, abs(O.DCS_shape_CALC), amplitude_threshold));
        set(O.plot, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        if test == 1
            % plotting reconstructed shape amplitude
            O.plot_true = patch(isosurface(O.N1grid*O.p_sam, O.N2grid*O.p_sam, O.N3grid*O.p_sam, abs(O.DCS_shape_REC), amplitude_threshold));
            set(O.plot_true, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
    else
        % plotting calculated shape isosurface
        [faces,verts,colors] = isosurface(O.N1grid*O.p_sam, O.N2grid*O.p_sam, O.N3grid*O.p_sam, abs(O.DCS_shape_CALC), amplitude_threshold, angle(O.DCS_shape_CALC));
        O.plot = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
        if test == 1
            % plotting reconstructed shape isosurface
            [faces,verts,colors] = isosurface(O.N1grid*O.p_sam, O.N2grid*O.p_sam, O.N3grid*O.p_sam, abs(O.DCS_shape_REC), amplitude_threshold, angle(O.DCS_shape_REC));
            O.plot_true = patch('Vertices', verts, 'Faces', faces, 'FaceVertexCData', colors, 'FaceColor', 'interp', 'edgecolor', 'none');
        end
    end
    
    % plotting x', y' and z' axes
    x_1p_axis = quiver3(0, 0, 0, 0.9*O.N1*O.p_sam/2, 0, 0);
    set(x_1p_axis, 'Color', 'black', 'Linewidth', 2, 'AutoScale', 'off');
    text(O.N1*O.p_sam/2, 0, 0, 'x''', 'Color', 'black', 'FontSize', 14);
    x_2p_axis = quiver3(0, 0, 0, 0, 0.9*O.N2*O.p_sam/2, 0);
    set(x_2p_axis, 'Color', 'black', 'Linewidth', 2, 'AutoScale', 'off');
    text(0, O.N2*O.p_sam/2, 0, 'y''', 'Color', 'black', 'FontSize', 14);
    x_3p_axis = quiver3(0, 0, 0, 0, 0, 0.9*O.N3*O.p_sam/2);
    set(x_3p_axis, 'Color', 'black', 'Linewidth', 2, 'AutoScale', 'off');
    text(0, 0, O.N3*O.p_sam/2, 'z''', 'Color', 'black', 'FontSize', 14);
    
    % overlap textbox
    if test == 1
        % centring masks for overlap calculation
        O.DCS_shape_CALC_MASK = circshift(O.DCS_shape_CALC_MASK, size(O.DCS_shape_CALC_MASK)/2-O.DCS_shape_CALC_COM);
        O.DCS_shape_REC_MASK = circshift(O.DCS_shape_REC_MASK, size(O.DCS_shape_REC_MASK)/2-O.DCS_shape_REC_COM);
        % calculating overlap
        overlap = round(abs((1-abs(sum(sum(sum(O.DCS_shape_CALC_MASK - O.DCS_shape_REC))))/sum(sum(sum(O.DCS_shape_REC))))*100), 2);
        annotation('textbox',[0.17, 0.1, .3, .3], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','String',['Overlap: ',num2str(overlap),'%'], 'BackgroundColor', 'white','FitBoxToText','on');
        legend([O.plot, O.plot_true], ['Calculated (', num2str(O.hkl(1)),num2str(O.hkl(2)),num2str(O.hkl(3)),') DCS shape'], ['Reconstructed (', num2str(O.hkl(1)),num2str(O.hkl(2)),num2str(O.hkl(3)),') DCS shape']);
    else
        legend([O.plot], ['Calculated (', num2str(O.hkl(1)),num2str(O.hkl(2)),num2str(O.hkl(3)),') DCS shape']);
    end
    
    % setting plot parameters
    title(['Calculated DCS shape vs. Reconstructed DCS shape for (', num2str(O.hkl(1)),num2str(O.hkl(2)),num2str(O.hkl(3)),')']);
    xlabel('x'' (m)'); ylabel('y'' (m)'); zlabel('z'' (m)');
    set(gca,'XDir','normal');
    set(gca,'YDir','normal');
    daspect([1,1,1]);
    axis equal;
    axis vis3d xy;
    grid on;
    xlim([-O.N1*O.p_sam/2, O.N1*O.p_sam/2]); ylim([-O.N2*O.p_sam/2, O.N2*O.p_sam/2]); zlim([-O.N3*O.p_sam/2, O.N3*O.p_sam/2]);
    view(viewpoint(1), viewpoint(2));
    lighting gouraud;
    if test ~= 1
        camlight('headlight');
    end
    if amplitudes ~= 1
        c = colorbar;
        ylabel(c, 'Phase');
        caxis([-pi, pi]);
    end
end

fprintf('\n...done\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = centerOfMass(A,varargin)
% CENTEROFMASS finds the center of mass of the N-dimensional input array
%
%   CENTEROFMASS(A) finds the gray-level-weighted center of mass of the
%   N-dimensional numerical array A. A must be real and finite. A warning
%   is issued if A contains any negative values. Any NaN elements of A will
%   automatically be ignored. CENTEROFMASS produces center of mass
%   coordinates in units of pixels. An empty array is returned if the
%   center of mass is undefined.
%
%   The center of mass is reported under the assumption that the first
%   pixel in each array dimension is centered at 1.
%
%   Also note that numerical arrays other than DOUBLE and SINGLE are
%   converted to SINGLE in order to prevent numerical roundoff error.
%
%   Examples:
%       A = rgb2gray(imread('saturn.png'));
%       C = centerOfMass(A);
%
%       figure; imagesc(A); colormap gray; axis image
%       hold on; plot(C(2),C(1),'rx')
%
%   See also: 
%
%

%
%   Jered R Wells
%   2013/05/07
%   jered [dot] wells [at] gmail [dot] com
%
%   v1.0
%
%   UPDATES
%       YYYY/MM/DD - jrw - v1.1
%
%

%% INPUT CHECK
narginchk(0,1);
nargoutchk(0,1);
fname = 'centerOfMass';

% Checked required inputs
validateattributes(A,{'numeric'},{'real','finite'},fname,'A',1);

%% INITIALIZE VARIABLES
A(isnan(A)) = 0;
if ~(strcmpi(class(A),'double') || strcmpi(class(A),'single'))
    A = single(A);
end
if any(A(:)<0)
    warning('MATLAB:centerOfMass:neg','Array A contains negative values.');
end

%% PROCESS
sz = size(A);
nd = ndims(A);
M = sum(A(:));
C = zeros(1,nd);
if M==0
    C = [];
else
    for ii = 1:nd
        shp = ones(1,nd);
        shp(ii) = sz(ii);
        rep = sz;
        rep(ii) = 1;
        ind = repmat(reshape(1:sz(ii),shp),rep);
        C(ii) = sum(ind(:).*A(:))./M;
    end
end

% Assemble the VARARGOUT cell array
varargout = {C};

end % MAIN