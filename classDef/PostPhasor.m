classdef PostPhasor < handle
    %   PHASOR Summary of this class goes here
    %   3D version
    % h - horizontal axis of the diffraction pattern
    % v - vertical axis of the diffraction pattern
    % r - rocking curve direction (3 dimension in arrays)
    % List of methods:
   
    properties
        data;
        dataTime; 
        derivatives;
        pre_path;
        path;
        experiment;                        
        object_fft_mod;    
        object_sampling;        
        object;
        displacement;
        strain;
        strainLab;
        strain_histogram;
        strain_histogram_vector;
        strain_mask;
        strain_mask_bulk;
        strain_mask_shell;
        support;
        volume;
        plotting;
        mask;           
        prtf;
        
        % Constants. They are needed for correct labeling of axes
        h                       = 4.1357e-15;                                  % Plank's constant
        c                       = 2.99792458e8;                                % Speed of light in vacuum
    end
        
    % - Constructor of the object -    
    methods
        function postPhasor = PostPhasor(input_param)           
            postPhasor.pre_path = input_param.pre_path;                       
            
            try
                postPhasor.object = input_param.object;
                postPhasor.mask = ones(size(postPhasor.object));
            catch 
                error('No reconstruction found!');
            end
            
            try
                postPhasor.experiment = input_param.experiment;
                postPhasor.experiment.wavelength = postPhasor.h*postPhasor.c/postPhasor.experiment.energy;
            catch
                warning('No experimental parameters found! Coordinate transform is not possible!')
            end
            
            try
                postPhasor.object_sampling = input_param.object_sampling;
                postPhasor.update_plotting_vectors;
            catch                
                error('No real-space sampling defined!');
            end
                
            try
                postPhasor.plotting.binning = input_param.binning;
                postPhasor.update_plotting_vectors;
            catch
                postPhasor.plotting.binning = [1,1];
                warning('No binning defined! Using [1,1]');
            end
            disp('PostPhasor instance created succesfully');
        end
    end  
    
    methods
        % Load methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function load_data_mat(postPhasor)
            % Dynamic determination of the dataset name
            in = load(postPhasor.experiment.data_path); 
            names = fieldnames(in);
            temp = in.(names{1});
            clear in;           
            
%             temp(isnan(temp)) = 0;
%             if temp<0
%                 postPhasor.data_gap = (temp==-1);
%                 temp(postPhasor.data_gap) = 0;
%                 disp('Gaps are recorded!');
%             else
%                 postPhasor.data_gap = (temp==-1);
%                 disp('No gaps');
%             end            
            
            postPhasor.data = single((double(temp)));
                        
            fprintf('Data loaded: %s\n',postPhasor.experiment.data_path);
            fprintf('Data size: [%d %d %d]\n',size(postPhasor.data));
        end    
        
        function load_object_mat(postPhasor)
            % Dynamic determination of the dataset name
            in = load(postPhasor.data_meta.reconstruction_path); 
            names = fieldnames(in);
            temp = in.(names{1});
            clear in;                                               
            
            postPhasor.object = temp;            
            
            fprintf('Sample reconstruction loaded: %s\n',postPhasor.data_meta.reconstruction_path);            
        end     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        % Class methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                         
        function update_plotting_vectors(postPhasor)
            postPhasor.plotting.object.vector1 = (-size(postPhasor.object,1)/2:size(postPhasor.object,1)/2-1).*postPhasor.object_sampling(1)*1e9; % [nm convention]
            postPhasor.plotting.object.vector2 = (-size(postPhasor.object,2)/2:size(postPhasor.object,2)/2-1).*postPhasor.object_sampling(2)*1e9;
            postPhasor.plotting.object.vector3 = (-size(postPhasor.object,3)/2:size(postPhasor.object,3)/2-1).*postPhasor.object_sampling(3)*1e9;
            
            if ~isempty(postPhasor.strain)
                postPhasor.plotting.strain.vector1 = (-size(postPhasor.strain,1)/2:size(postPhasor.strain,1)/2-1).*postPhasor.object_sampling(1)*1e9; % [nm convention]
                postPhasor.plotting.strain.vector2 = (-size(postPhasor.strain,2)/2:size(postPhasor.strain,2)/2-1).*postPhasor.object_sampling(2)*1e9;
                postPhasor.plotting.strain.vector3 = (-size(postPhasor.strain,3)/2:size(postPhasor.strain,3)/2-1).*postPhasor.object_sampling(3)*1e9;
                disp('Plotting vectors for strain updated!');
            end
            
            disp('Plotting vectors for object updated!');
        end
        
        function create_mask(postPhasor,value)    
            if nargin == 1
                value = 0.1;
            end
            postPhasor.mask = abs(postPhasor.object) > value*max(abs(postPhasor.object(:)));
            fprintf('Mask is updated to %.2f of max amplitude\n',value);
        end
        
        function calculate_object_fft_mod(postPhasor)
            % Should we compare the same volume inside the support
                % only, since everything outside is not an object?
            postPhasor.object_fft_mod = abs(fftNc(postPhasor.object));            
        end                                                            
        
        % Post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function flip(postPhasor,dimension)
            try
                if nargin == 2
                    postPhasor.object = flip(postPhasor.object,dimension);
                    postPhasor.mask = flip(postPhasor.mask,dimension);
                    fprintf('Object is flipped along dimension %d\n',dimension);
                    if ~isempty(postPhasor.strain)
                        postPhasor.strain = flip(postPhasor.strain);
                        postPhasor.strain_mask = flip(postPhasor.strain_mask);
                        fprintf('Strain is flipped along dimension %d\n',dimension);
                    end
                    if ~isempty(postPhasor.displacement)
                        postPhasor.displacement = flip(postPhasor.strain);
                        postPhasor.strain_mask = flip(postPhasor.strain_mask);
                        fprintf('Strain is flipped along dimension %d\n',dimension);
                    end
                else
                    postPhasor.object = flip(postPhasor.object);
                    postPhasor.mask = flip(postPhasor.mask);
                    fprintf('Object is flipped along all dimensions\n');
                    if ~isempty(postPhasor.strain)
                        postPhasor.strain = flip(postPhasor.strain);
                        postPhasor.strain_mask = flip(postPhasor.strain_mask);
                        fprintf('Strain is flipped along all dimensions');
                    end 
                end
            catch
                error('Object can not be flipped along selected dimension');
            end
            
        end             
        
        function crop_manual(postPhasor)  
            
            figure; 
            imagesc(squeeze(sum(abs(postPhasor.object),3)));
            hIm = imrect;
            pos = round(getPosition(hIm)); %[xmin ymin width height]                        
            
            maskFlag = 0;
            strainFlag = 0;
            dispFlag = 0;
            
            for ii = 1:size(postPhasor.object,3)      
                temp.object(:,:,ii) = squeeze(postPhasor.object(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),ii)); 
                try
                    temp.mask(:,:,ii) = squeeze(postPhasor.mask(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),ii));
                    maskFlag = 1;
                catch
                    maskFlag = 0;
                end

                try
                    temp.strain(:,:,ii) = squeeze(postPhasor.strain(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),ii));
                    strainFlag = 1;
                catch
                    strainFlag = 0;
                end

                try
                    temp.displacement(:,:,ii) = squeeze(postPhasor.displacement(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),ii));
                    dispFlag = 1;
                catch
                    dispFlag = 0;
                end
            end
            
            postPhasor.object = temp.object;
            
            if maskFlag            
                postPhasor.mask = temp.mask;
            end
            
            if strainFlag
                postPhasor.strain = temp.strain;
            end
            
            if dispFlag
                postPhasor.displacement = temp.displacement;
            end
            
            clear temp
            close;
           
            figure; 
            imagesc(squeeze(sum(abs(postPhasor.object),2)));
            hIm = imrect;
            pos = round(getPosition(hIm)); %[xmin ymin width height]                        

            for ii = 1:size(postPhasor.object,2)      
                temp.object(:,ii,:) = squeeze(postPhasor.object(pos(2):pos(2)+pos(4),ii,pos(1):pos(1)+pos(3)));
                try
                    temp.mask(:,ii,:) = squeeze(postPhasor.mask(pos(2):pos(2)+pos(4),ii,pos(1):pos(1)+pos(3)));   
                catch
                end
                
                try
                    temp.strain(:,ii,:) = squeeze(postPhasor.strain(pos(2):pos(2)+pos(4),ii,pos(1):pos(1)+pos(3)));
                catch                    
                end

                try
                    temp.displacement(:,ii,:) = squeeze(postPhasor.displacement(pos(2):pos(2)+pos(4),ii,pos(1):pos(1)+pos(3)));
                catch                    
                end
            end
            
            postPhasor.object = temp.object;
            if maskFlag            
                postPhasor.mask = temp.mask;
            end
            
            if strainFlag
                postPhasor.strain = temp.strain;
            end
            
            if dispFlag
                postPhasor.displacement = temp.displacement;
            end
            
            clear temp
            close;
            
            postPhasor.update_plotting_vectors;
            
        end
        
        function apodize(postPhasor,window_type,window_width)
            % In frequency domain
%             if nargin == 1
%                 window_type = 'Tukey';
%                 window_width = 0.9;            
%             end
            
            switch window_type
                case 'Tukey'
                    tukeyx=tukeywin(size(postPhasor.object,1),window_width); %1D window
                    tukeyy=tukeywin(size(postPhasor.object,2),window_width);
                    tukeyz=tukeywin(size(postPhasor.object,3),window_width);

                    tukey2=ones(size(postPhasor.object,1),size(postPhasor.object,2));
                    window=ones(size(postPhasor.object,1),size(postPhasor.object,2),size(postPhasor.object,3));

                    for i=1:size(postPhasor.object,1)              % there is probably a much better way to compute this
                        tukey2(i, :)=tukeyx(i).*tukeyy;
                        for j=1:size(postPhasor.object,2)
                            window(i, j, :)=tukey2(i, j).*tukeyz;
                        end
                    end
            end
            fft_object = fftNc(postPhasor.object);
            fft_max = max(abs(fft_object(:)));
            fft_window = fftNc(postPhasor.object).*window;
            fft_window = fft_window.*fft_max./max(abs(fft_window(:)));
            postPhasor.object = ifftNc(fft_window);   
            fprintf('Object is apodized by %s window of %.2f width\n',window_type,window_width);
        end
        
        function center_data(postPhasor)
            if ~isempty(postPhasor.data)
                [postPhasor.data,center] = center_array_com(postPhasor.data);
                fprintf('Data centered at [%.2f %.2f %.2f]\n',(center));
            else
                error('No data found!')
            end
        end                
        
        function matchSampling(postPhasor)
            if postPhasor.object_sampling(1)~=postPhasor.object_sampling(2)~=postPhasor.object_sampling(3)
                min_sampling = min(postPhasor.object_sampling);
                
                ratio = postPhasor.object_sampling/min_sampling;
                                
                if ratio(1) == 1
                    for ii = 1:size(postPhasor.object,1)
                        t(ii,:,:) = imresize(squeeze(postPhasor.object(ii,:,:)),[size(postPhasor.object,2)*ratio(2),size(postPhasor.object,3)*ratio(3)]);
                    end
                elseif ratio(2) == 1
                    for ii = 1:size(postPhasor.object,2)
                        t(:,ii,:) = imresize(squeeze(postPhasor.object(:,ii,:)),[size(postPhasor.object,1)*ratio(1),size(postPhasor.object,3)*ratio(3)]);
                    end
                elseif ratio(3) == 1
                    for ii = 1:size(postPhasor.object,3)
                        t(:,:,ii) = imresize(squeeze(postPhasor.object(:,:,ii)),[size(postPhasor.object,1)*ratio(1),size(postPhasor.object,2)*ratio(2)]);
                    end
                end
                                
                postPhasor.object = t;
                postPhasor.object_sampling = [min_sampling,min_sampling,min_sampling];
                postPhasor.update_plotting_vectors;
                disp('+ Sampling matched!')
            else
                disp('Sampling is already matched!')
            end
        end
        
        function center_object(postPhasor)
            [postPhasor.object, com_val] = center_array_com(postPhasor.object);            
            fprintf('+ Object centered at [%.1f %.1f %.1f]!\n',com_val(1),com_val(2),com_val(3));
        end
        
        function rotateObject(postPhasor,rotationAngles)
            try
                [x, y, z] = meshgrid(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,postPhasor.plotting.object.vector3);                       

                U = rotxd(rotationAngles(1))*rotyd(rotationAngles(2))*rotzd(rotationAngles(3));

                coordinatesNew = [x(:),y(:),z(:)]*U;
                postPhasor.object = interp3(x, y, z, postPhasor.object, reshape(coordinatesNew(:,1),size(postPhasor.object)), reshape(coordinatesNew(:,2),size(postPhasor.object)),reshape(coordinatesNew(:,3),size(postPhasor.object)), 'linear', 0); % make any values outside original data zero.
%                 postPhasor.object_sampling = [abs(coordinatesNew(2,1)-coordinatesNew(1,1)),abs(coordinatesNew(2,2)-coordinatesNew(1,2)),abs(coordinatesNew(2,3)-coordinatesNew(1,3))];                
                postPhasor.update_plotting_vectors;
                disp('+ Object rotated!')
            catch
                warning('Can not rotate object!');
            end
        end
            
        function remove_phase_ramp(postPhasor)
            c0 = size(postPhasor.object)/2;                                                   
            if numel(size(postPhasor.object)) == 3 
                for ii = 1:3
                    dummy = fftNc(postPhasor.object);
                    c1 = ndimCOM(abs(dummy).^4,'auto');                
                    d =((-1)^ii)*(c1-c0);                                       
                    fprintf('Data center: [%.2f, %.2f, %.2f]\nShift by [%.2f, %.2f, %.2f]\n', c1(1), c1(2), c1(3),d(1),d(2),d(3));                    
                    dummy = imtranslate(dummy,d);
                    dummy = ifftNc(dummy);
                    % Average phase value
                    postPhasor.object = dummy.*exp(-1j*mean(angle(dummy(:))));
                end
            elseif numel(postPhasor.data_meta.data_size) == 2
                c1 = ndimCOM(postPhasor.data,'auto');
                fprintf('Data center: [%.2f, %.2f]', c1(1), c1(2));
            else
                warning('Data dimensionality is wrong!');
            end                    
            disp('Phase ramp removed!');
        end
        
%         function transform2lab_backup(postPhasor)
%             % Unified coordinate transformation
%             DCS_to_SS_backup;                             
%             postPhasor.update_plotting_vectors;
%         end
        
        function transform2lab(postPhasor)
            %change sign of the phase to have a correct displacement field
%             postPhasor.object = abs(postPhasor.object).*exp(-1j.*angle(postPhasor.object));
            % Unified coordinate transformation           
            
            DCS_to_SS_LAB; 
            postPhasor.update_plotting_vectors;
        end
        
        function transform2sample(postPhasor)
            %change sign of the phase to have a correct displacement field
%             postPhasor.object = abs(postPhasor.object).*exp(-1j.*angle(postPhasor.object));
            % Unified coordinate transformation           
            
            DCS_to_SS; 
            postPhasor.update_plotting_vectors;
        end
        
        function add_pi(postPhasor)
            postPhasor.object = postPhasor.object.*exp(1j*pi);
            disp('+pi value is added to the phase of the object')
        end
        
        function add_half_pi(postPhasor)
            postPhasor.object = postPhasor.object.*exp(1j*pi/2);
            disp('+pi/2 value is added to the phase of the object')
        end
        
        function twin_object(postPhasor)
            F = ifftshift(fftn(fftshift(postPhasor.object)));
            postPhasor.object = fftshift(ifftn(ifftshift(conj(F))));
        end                
        
        function calculate_qVectorLab(postPhasor)
            % For beamline geometries go to: https://github.com/DzhigaevD/phasor/wiki/Beamline-geometries
            % Modulus of the incidence wave-vector
            k_mod = (2*pi)/postPhasor.experiment.wavelength;
            
            % The direction of the incidence beam in the lab frame
            postPhasor.experiment.k_i = [0;0;k_mod];
            
            try
                switch postPhasor.experiment.beamline
                    % Unified coordinate transformation
                    case '34idc'  
                        % delta - positive rotation around vertical axis
                        % gamma - negative rotation around horizontal axis

                        % The direction of the scattered beam in the lab frame
                        % the + sign before the sind (2nd element) of the
                        % matrix is there due to the rotation direction of
                        % the detector arm                                              
                        postPhasor.experiment.k_s = [sind(postPhasor.experiment.delta+postPhasor.experiment.delta_correction)*cosd(postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction);...
                        sind(postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction);...
                        cosd(postPhasor.experiment.delta+postPhasor.experiment.delta_correction)*cosd(postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction)].*k_mod;
                        
                    case 'p10'
                        % gamma - positive rotation around vertical axis
                        % delta - negative rotation around horizontal axis
                        postPhasor.experiment.k_s = [sind(postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction)*cosd(postPhasor.experiment.delta+postPhasor.experiment.delta_correction);...
                        sind(postPhasor.experiment.delta+postPhasor.experiment.delta_correction);...
                        cosd(postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction)*cosd(postPhasor.experiment.delta+postPhasor.experiment.delta_correction)].*k_mod;
                    
                    case 'nanomax'
                        % gamma - negative rotation around vertical axis
                        % delta - negative rotation around horizontal axis
                        postPhasor.experiment.k_s = [sind(-postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction)*cosd(postPhasor.experiment.delta+postPhasor.experiment.delta_correction);...
                        sind(postPhasor.experiment.delta+postPhasor.experiment.delta_correction);...
                        cosd(-postPhasor.experiment.gamma+postPhasor.experiment.gamma_correction)*cosd(postPhasor.experiment.delta+postPhasor.experiment.delta_correction)].*k_mod;                          
                end
                
                % The direction of the wave transfer in the lab frame
                postPhasor.experiment.qVectorLab = postPhasor.experiment.k_s-postPhasor.experiment.k_i;
                fprintf('+ Calculated Q-vector in lab system of %s!\n',postPhasor.experiment.beamline);
                
                disp('Derivatives from the geometry:');
                twoTheta = acosd(dot(postPhasor.experiment.k_i,postPhasor.experiment.k_s)/(norm(postPhasor.experiment.k_i)*norm(postPhasor.experiment.k_s)));
                fprintf('2Theta is: %.3f\n',twoTheta);
                fprintf('|Q| is: %.3f\n',norm(postPhasor.experiment.qVectorLab)*1e-10);
                fprintf('d-spacing is: %.3f\n',2*pi/norm(postPhasor.experiment.qVectorLab)*1e10);
            catch
                error('Can not calculate Q-vector in lab system for %s!',postPhasor.experiment.beamline);
            end           
        end
        
        function calculate_volume(postPhasor)
            if ~isempty(postPhasor.support) 
                postPhasor.volume = sum(postPhasor.support(:))*prod(postPhasor.object_sampling.*1e9);
                fprintf('The volume calculated from support: %f nm³\n', postPhasor.volume);
            else
                threshold = 0.3;           
                v = abs(postPhasor.object) > threshold;
                postPhasor.volume = sum(v(:))*prod(postPhasor.object_sampling.*1e9);
                fprintf('The volume calculated from %.1f of amplitude: %f nm³\n', threshold, postPhasor.volume);
            end                        
        end
        
        function calculate_displacement(postPhasor)
        % displacment field: phase devided by modulus q
            try
                postPhasor.displacement = -angle(postPhasor.object)/norm(postPhasor.experiment.qVectorLab); % displacment field
                disp('+ Displacement field is calculated as negative phase, for MATLAB FFT')
            catch
                try
                    postPhasor.calculate_qVectorLab;
                    postPhasor.calculate_displacement;
                catch
                    error('Can not calculate displacement, check experimental parameters!');
                end
            end
        end                
        
        function calculate_strainLab(postPhasor)
            try
                % strain calc   

                % vector in the q deriction; length in inv A
                qplot = [postPhasor.experiment.qVectorLab(1); postPhasor.experiment.qVectorLab(2); postPhasor.experiment.qVectorLab(3)];

                % unit vector in the direction of q
                qplot_unitary = qplot./sqrt(qplot(1).^2+qplot(2).^2+qplot(3).^2);%      

                % displacment field: phase devided by modulus q
                postPhasor.calculate_displacement;

                % gradient of phase in x y z with spacing from interpolation
                [e_x, e_y, e_z] = gradient(postPhasor.displacement, postPhasor.object_sampling(1),postPhasor.object_sampling(2),postPhasor.object_sampling(3)); %e_x: du/dx
                
                % sum of the progection of gradient on the q direction                
                postPhasor.strain = (e_x*qplot_unitary(1) + e_y*qplot_unitary(2) + e_z*qplot_unitary(3));
                postPhasor.strain_mask = postPhasor.mask;
                postPhasor.update_plotting_vectors;
            catch
                try
                    postPhasor.calculate_qVectorLab;
                    postPhasor.calculate_strainLab;
                catch
                    error('Can not calculate strain, check experimental parameters!');
                end
            end
        end
        
%         function calculate_strain(postPhasor, strain_axis)
%             warning('This is a depricated version of strain calculation, will be removed soon! Use calculate_strainLab!');
%             H = 2*pi/postPhasor.experiment.d_spacing;
% 
%             postPhasor.displacement = angle(postPhasor.object)./H;
%             
%             if strain_axis < 0
%                 postPhasor.strain = flip(diff(flip(postPhasor.displacement,abs(strain_axis)),1,abs(strain_axis))./postPhasor.object_sampling(abs(strain_axis)),abs(strain_axis));
%             else
%                 postPhasor.strain = diff(postPhasor.displacement,1,strain_axis)./postPhasor.object_sampling(strain_axis);
%             end
%             
%             mask_shift = [0,0,0];
%             mask_shift(abs(strain_axis)) = 1;
% 
%             if abs(strain_axis) == 1
%                 postPhasor.strain_mask = postPhasor.mask(1:end-1,:,:);
%             elseif abs(strain_axis) == 2
%                 postPhasor.strain_mask = postPhasor.mask(:,1:end-1,:);
%             elseif abs(strain_axis) == 3
%                 postPhasor.strain_mask = postPhasor.mask(:,:,1:end-1);
%             end
% 
%             postPhasor.strain_mask = postPhasor.strain_mask+circshift(postPhasor.strain_mask,mask_shift);
%             postPhasor.strain_mask = postPhasor.strain_mask == 2;
%             
%             postPhasor.update_plotting_vectors;
%         end       
            
        function calculate_strain_unwrap(postPhasor, strain_axis)
            %jclark
            %unwrap based on diff and cumsum
            %drc is direction

            if exist('strain_axis') ~= 1
                strain_axis = 1;
            end
            strain_axis_unwrap = abs(strain_axis);
            %pad the dimension
            temp=exp(1j*pad_by_one(angle(postPhasor.object),strain_axis_unwrap));

            switch strain_axis_unwrap
                case 1
                    phdiff=circshift(temp,-1).*conj(temp);
                case 2
                    phdiff=circshift(temp,[0 -1]).*conj(temp);
                case 3
                    phdiff=circshift(temp,[0 0 -1]).*conj(temp);
                case 4
                    phdiff=circshift(temp,[0 0 0 -1]).*conj(temp);
            end

            %calc the derivative of phase
            temp=angle(phdiff);

            %now do the cumsum
            temp=cumsum(temp,strain_axis_unwrap);
            
            switch strain_axis_unwrap
                case 1
                    temp=circshift(temp,1);
                case 2

                    temp=circshift(temp,[0,1]);
                case 3

                    temp=circshift(temp,[0,0,1]);
                case 4

                    temp=circshift(temp,[0,0,0,1]);
            end
            %now get original size
            phase = crop_by_one(temp,strain_axis_unwrap);
            
            H = 2*pi/postPhasor.experiment.d_spacing;
            postPhasor.displacement = -phase./H;
            
            if strain_axis < 0
                postPhasor.strain = flip(diff(flip(postPhasor.displacement,abs(strain_axis)),abs(strain_axis))./postPhasor.object_sampling(strain_axis));
            else
                postPhasor.strain = flip(diff(postPhasor.displacement,strain_axis)./postPhasor.object_sampling(strain_axis));
            end
            
            mask_shift = [0,0,0];
            mask_shift(abs(strain_axis)) = 1;

            if abs(strain_axis) == 1
                postPhasor.strain_mask = postPhasor.mask(1:end-1,:,:);
            end

            postPhasor.strain_mask = postPhasor.strain_mask+circshift(postPhasor.strain_mask,mask_shift);
            postPhasor.strain_mask = postPhasor.strain_mask == 2;
            
            postPhasor.update_plotting_vectors;
        end
        
        function calculate_strain_histogram(postPhasor, domain)
            % Calculation of histogram inside the domain selected
            % Plot the histogram
            
            if ~exist('domain')
                domain = 'full';                
            end
            % Use an input parameter to show other complex valued matrix
            switch domain
                case 'full'
                    strain_mask = postPhasor.strain_mask;
                    title_s = 'Full object';
                    disp('Histogram of strain in the full object volume');
                case 'bulk'
                    if isempty(postPhasor.strain_mask_bulk)
                        error('No segmentation done! Use .segment_strain_mask')
                    else
                        strain_mask = postPhasor.strain_mask_bulk;
                        title_s = 'Bulk';
                        disp('Histogram of strain in the bulk object volume');                    
                    end
                case 'shell'
                    if isempty(postPhasor.strain_mask_shell)
                        error('No segmentation done! Use .segment_strain_mask')
                    else
                        strain_mask = postPhasor.strain_mask_shell;
                        title_s = 'Shell';
                        disp('Histogram of strain in the shell object volume');
                    end
                case 'both'
                    if isempty(postPhasor.strain_mask_bulk)
                        error('No segmentation done! Use .segment_strain_mask')
                    else
                        strain_mask(:,:,:,1) = postPhasor.strain_mask_bulk;
                        strain_mask(:,:,:,2) = postPhasor.strain_mask_shell;                        
                        title_s = 'Bulk and shell';
                        legend_s{1} = 'Bulk';
                        legend_s{2} = 'Shell';
                        disp('Histogram of strain in the bulk and shell of the object volume');                    
                    end
            end
            
            for ii = 1:size(strain_mask,4)
                input(:,:,:,ii) = postPhasor.strain.*strain_mask(:,:,:,ii);
            end            
                        
            figure;
            for ii = 1:size(strain_mask,4)
                temp = input(:,:,:,ii);
                
                [hV,edges] = histcounts(temp,200,'Normalization','probability'); % 
                max_val = max(hV);
                [val, pos] = find(hV == max_val);
                fprintf('%.2f%% of strain values are in the range: [%.2e : %.2e]%%\n', max(hV)*100, edges(pos)*100, edges(pos+1)*100);                      
                                
                hH = histogram(temp(temp~=0),1000,'Normalization','probability','EdgeColor','none');  hold on;                
            end
            
            set(gca,'FontSize',24); 
            set(gcf,'Units','normalized','Position',[0.1,0.1,0.5,0.5]); 
            xline(0,'--');
            yline(max(hH.Values(:))/2,'r'); 
            yline(max(hH.Values(:))/4,'g');             
            xlabel('Strain');
            ylabel('Probability');         
            title(title_s);
            
            if ii == 2
                legend(legend_s{1},legend_s{2},'Zero strain','FWHM','2FWHM')
            else
                legend(title_s,'Zero strain','FWHM','2FWHM')
            end
            
%             postPhasor.strain_histogram = hH.Values;
%             postPhasor.strain_histogram_vector = hH.BinEdges(1:end-1);                        
            
            
%             figure; 
%             plot(hH.BinEdges(1:end-1),log10(hH.Values));
%             title('log-plot of probability');
        end

        function calculate_central_profile_amplitude(postPhasor)                                       
            val = squeeze(abs(postPhasor.object(:,round(end/2),round(end/2))));
            val_contour = squeeze((postPhasor.mask(:,round(end/2),round(end/2)))).*max(val(:));
            val_contour_mask = val_contour>0;
            val_contour_length1 = sum(val_contour_mask).*postPhasor.object_sampling(1)*1e9;
            postPhasor.derivatives.central_profile_amplitude1 = val;

            val = squeeze(abs(postPhasor.object(round(end/2),:,round(end/2))));
            val_contour = squeeze((postPhasor.mask(round(end/2),:,round(end/2)))).*max(val(:));                
            val_contour_mask = val_contour>0;
            val_contour_length2 = sum(val_contour_mask).*postPhasor.object_sampling(2)*1e9;
            postPhasor.derivatives.central_profile_amplitude2 = val;

            val = squeeze(abs(postPhasor.object(round(end/2),round(end/2),:)));
            val_contour = squeeze(postPhasor.mask(round(end/2),round(end/2),:)).*max(val(:));
            val_contour_mask = val_contour>0;
            val_contour_length3 = sum(val_contour_mask).*postPhasor.object_sampling(3)*1e9;
            postPhasor.derivatives.central_profile_amplitude3 = val;

            postPhasor.derivatives.central_profile_lengths = [val_contour_length1, val_contour_length2, val_contour_length3];
            
            disp('+ Object shape profiles are saved to derivatives in postPhasor!');
        end
        
        function calculate_central_profile_strain(postPhasor)     
            try
                val = squeeze(postPhasor.strain(:,round(end/2),round(end/2)));
                val_contour = squeeze((postPhasor.mask(:,round(end/2),round(end/2)))).*max(val(:));
                val_contour_mask = val_contour>0;
                postPhasor.derivatives.central_profiles_strain1 = val.*val_contour_mask;

                val = squeeze(abs(postPhasor.strain(round(end/2),:,round(end/2))));
                val_contour = squeeze((postPhasor.mask(round(end/2),:,round(end/2)))).*max(val(:));                
                val_contour_mask = val_contour>0;
                postPhasor.derivatives.central_profiles_strain2 = val.*val_contour_mask;
                
                val = squeeze(abs(postPhasor.strain(round(end/2),round(end/2),:)));
                val_contour = squeeze(postPhasor.mask(round(end/2),round(end/2),:)).*max(val(:));
                val_contour_mask = val_contour>0;
                postPhasor.derivatives.central_profiles_strain3 = val.*val_contour_mask;       
                
                disp('+ Object strain profiles are saved to derivatives in postPhasor!');
            catch
                postPhasor.calculate_qVectorLab;
                postPhasor.calculate_strainLab;
                calculate_central_profile_strain;
            end
        end
        
        function calculate_angular_profile(postPhasor, source, slice, center, radius, startAngle)
            startAngle = startAngle*pi/180;
            switch source
                case 'displacement'
                    if slice(1)~=0
                        tImg = squeeze(postPhasor.displacement(slice(1),:,:));
                        tAmp = squeeze(abs(postPhasor.object(slice(1),:,:)));
                        axisN = 1;
                    elseif slice(2)~=0
                        tImg = squeeze(postPhasor.displacement(:,slice(2),:));
                        tAmp = squeeze(abs(postPhasor.object(:,slice(2),:)));                                                
                        axisN = 2;
                    elseif slice(3)~=0
                        tImg = squeeze(postPhasor.displacement(:,:,slice(3)));
                        tAmp = squeeze(abs(postPhasor.object(:,:,slice(3))));
                        axisN = 3;
                    end
                    labelY = 'Displacement,[m]';
                case 'phase'
                    tImg = squeeze(angle(postPhasor.object(slice(1),:,:)));  
                    labelY = 'Phase,[rad]';
            end
            figure;imagesc(tImg);axis image;colormap jet;
            
            [hd,vd] = meshgrid(-size(tImg,2)/2+1:(size(tImg,2)/2),-size(tImg,1)/2+1:(size(tImg,1)/2));

            hd = (hd+size(tImg,2)/2-center(1));
            vd = (vd+size(tImg,1)/2-center(2));
            
            t = 0:2*pi/18:2*pi;
            t = t+startAngle;
            xq = cos(t).*radius;
            yq = sin(t).*radius;
            
            val = interp2(hd,vd,tImg,xq',yq');
            
            handle = figure;
            ax1 = subplot(2,2,1); imagesc(hd(1,:),vd(:,1),tAmp);axis image; colorbar;
                title('Amplitude of the object');
            subplot(2,2,2); ax2 = imagesc(hd(1,:),vd(:,1),tImg);axis image; colormap jet; colorbar;
                title('Profile is calculated clockwise from dot marker')
                xline(0);
                yline(0);
                line(xq,yq,'Color',[1,1,1]);hold on;
                scatter(xq(1),yq(1),30,[0.25,0.25,0.25],'filled');
            ax3 = subplot(2,2,[3,4]); plot(t,val);
            xlabel('Angle, [rad]'); ylabel(labelY)     
            
            % GUI
%             uicontrol('Parent',handle,'Style',...
%             'slider','Min',1,'Max',size(postPhasor.object,axis),...
%             'Value',0,'Units','Normalized',...
%             'Position', [0.1 0.05 0.3 0.03],...
%             'Callback', @slideSlice);
        
%             function slideSlice(hObj,callbackdata)
%                 posVal = get(hObj,'Value');                 
%                 plotSlice(input,isoVal);
%             end
%             
%             function plotSlice(input, posVal)
%                 end
        end
        
        function segment_strain_mask(postPhasor,sigma,threshold)
            if nargin == 1
                sigma = 3;
                threshold = 0.75;
            elseif nargin == 2
                threshold = 0.75;
            end
            
            if isempty(postPhasor.object)
                error('No object defined! Create object first!')
            else
                try
                    temp = imgaussfilt3(double(postPhasor.strain_mask),[sigma sigma sigma]);
                catch
                    error('Wrong dimensionality of the input data. Expect 3D');
                end
                                
                postPhasor.strain_mask_bulk = temp>threshold; 
                postPhasor.strain_mask_shell = postPhasor.strain_mask-postPhasor.strain_mask_bulk; 
                          
                % Estimate shell thickness
                a = zeros(100,1);
                a(25:75) = 1;
                b = imgaussfilt(a,sigma);
                c = b>threshold;                
                s = sum(a-c)/2;
                th = s*postPhasor.object_sampling(1);
                fprintf('Mask is shrinked by sigma %.3f at threshold %.3f\n',sigma, threshold);
                fprintf('Shell thickness: %.2f nm\n',th*1e9);
            end
        end
        
        function calculate_prtf(postPhasor)
            % [IN DEVELOPMENT]
            if exist('np') ~= 1
                np=1;
            end
            
            if ~exist(postPhasor.experiment.data_path)
                error('No data path found!')
            else
                try
                    load_data_mat(postPhasor);
                catch
                    error('Can not load the data!')
                end
            end                           

            %get pixel sizes direction
            angxy=(postPhasor.experiment.sample_detector_d.*postPhasor.experiment.wavelength)./[size(postPhasor.data,1),size(postPhasor.data,2)]./postPhasor.experiment.detector_pitch;

            angz=(postPhasor.experiment.wavelength)./sind(size(postPhasor.data,3).*postPhasor.experiment.angular_step);

            rescaledxy=round(angxy/angz.*[size(postPhasor.data,1),size(postPhasor.data,2)]);

            FOVxy=angxy.*[size(postPhasor.data,1),size(postPhasor.data,2)];
            FOVz=angz*size(postPhasor.data,3);

            scalef=FOVxy/FOVz;

            newz=round(scalef*size(postPhasor.data,3));           

            Ipnm=(abs(fftshift(fftn(fftshift(postPhasor.object)))).^2);
%             Ipnm=convn(abs(fftshift(fftn(fftshift(pn)))).^2,params.coh,'same');

            Ipnm(Ipnm <0)=0;

            %get a threshold for fter shifting/sampling
            th=min(postPhasor.data(postPhasor.data ~= 0));

            %pp=0;
            %for qq=1:2:size(Ipnm,3)
            %    pp=pp+1;
            %    Ipnmn(:,:,pp)=sum(Ipnm(:,:,qq:qq+1),3);
            %    datan(:,:,pp)=sum(data(:,:,qq:qq+1),3);
            %end

            %fft and crop/pad to get same sampling
            %Ipnmn=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(Ipnm))),size(Ipnm,2),size(Ipnm,1),62)))));
            %datan=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(data))),size(Ipnm,2),size(Ipnm,1),62)))));

            xyz=center_of_mass(postPhasor.data);
            postPhasor.data=circshift(postPhasor.data,-round([xyz(2),xyz(1),xyz(3)]));
            Ipnm=circshift(Ipnm,-round([xyz(2),xyz(1),xyz(3)]));

            %change th FOV to be the same so the sampling is the same
            %Ipnmn=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(Ipnm))),size(Ipnm,2),size(Ipnm,1),newz(1))))));
            %datan=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(data))),size(Ipnm,2),size(Ipnm,1),newz(1))))));
            %now sampling is the same as dx dy
            %Ipnm=Ipnmn;
            %data=datan;

            %normalize
            Ipnm=Ipnm./sum(Ipnm(:))*sum(postPhasor.data(:));

            %sqrt of data
            Ipnm=sqrt(Ipnm);
            postPhasor.data=sqrt(postPhasor.data);

            %threshold
            postPhasor.data(postPhasor.data < th)=0;

            ratio=Ipnm./postPhasor.data;

            ratio(isnan(ratio)) = 0;
            ratio(isinf(ratio)) = 0;
            ratio(ratio > 2) = 1;

            ratio=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(ratio))),size(Ipnm,2),size(Ipnm,1),newz(1))))));

            nd=max(size(ratio));
            ratio=zero_pad_ver3(ratio,nd,nd,nd);

            ratio(ratio < .01)=0;
            ratio(isnan(ratio)) = 0;

            rad=1:(round(nd(1)));

            parfor qq=2:numel(rad)
                fprintf('Calculating shell %d out of %d\n',qq,numel(rad));
                [ circ ] = generate_sphere(nd(1),rad(qq))-generate_sphere(nd(1),rad(qq-1));

                circ0=circ.*(ratio > 0);

                ind=find(circ0 == 1);

                prtf(qq)=sum((ratio(ind)))./numel(ind(:));

%                 ind0=find(circ == 1);
%                 prt0(qq)=sum((ratio(ind0)))./numel(ind0(:));

            end
            postPhasor.prtf.values = prtf;
            postPhasor.prtf.q_spacing = 2*pi*postPhasor.experiment.detector_pitch./postPhasor.experiment.wavelength./postPhasor.experiment.sample_detector_d;
        end
               
        function calculate_resolution(postPhasor)
            % [IN DEVELOPMENT]
        end
        
        % Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_prtf(postPhasor)
            figure;
            try                
                plot(postPhasor.prtf);hold on;
                yline(1/exp(1));
            catch
                error('No prtf found!');
            end
        end
        
        function slice3d_object(postPhasor)
            vis3d(double(angle(postPhasor.object)), abs(postPhasor.object)>0.1);
        end
        
        function slice3dStrain(postPhasor)
            vis3d(double(postPhasor.strain), postPhasor.strain_mask);
        end
        
        function iso3dStrainDistribution(postPhasor,domain)
            if ~exist('domain')
                domain = 'full';
            end
            % Use an input parameter to show other complex valued matrix
            switch domain
                case 'full'
                    strain_mask = postPhasor.strain_mask;  
                    title_s = 'Full object';
                    disp('Displaying full strain isosurface');
                case 'bulk'
                    strain_mask = postPhasor.strain_mask_bulk; 
                    title_s = 'Bulk';
                    disp('Displaying bulk strain isosurface');
                case 'shell'
                    strain_mask = postPhasor.strain_mask_shell;  
                    title_s = 'Shell';
                    disp('Displaying shell strain isosurface');
            end
            input = postPhasor.strain.*strain_mask;            
            handle = figure;     
            set(gcf,'Units','Normalized','Position',[.1 0.1 .5 .8]);            
            panel = uipanel('Parent',handle,'Title','Strain distribution','FontSize',...
            12,'Units','Normalized','Position',[.1 0.1 .77 .85]); 
            
            ax = axes('Parent',panel); 
            title(title_s);
            
            uicontrol('Parent',handle,'Style',...
            'slider','Min',min(input(:)),'Max',max(input(:)),...
            'Value',0,'Units','Normalized',...
            'Position', [0.1 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReal); 
            
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.2,'Units','Normalized',...
            'Position', [0.7 0.05 0.2 0.03],...
            'Callback', @slideAlpha);
        
            isoVal = min(input(:)); 
            alphaVal = 0.2;
            hText = uicontrol('Parent',handle,'Style','text','String',sprintf('Strain value: %.4f',isoVal),'Units','Normalized',...
            'Position', [0.4 0.05 0.3 0.03]);
            drawIsosurface(input,alphaVal,isoVal,strain_mask);            
            
            function drawIsosurface(input,alphaVal,isoVal,strain_mask)
                cla(ax);
                axes(ax);
                
                % Shape isosurface  
                isosurface(strain_mask); alpha(alphaVal)
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
                hold on;            
                
                % Strain isosurface
                isosurface(input,isoVal);hold on                
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
                rotate3d on;
                grid on;
                axis tight;
                axis equal; 
                axis vis3d;
                axis image;
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];           
                colormap jet;    
                set(hText,'String',sprintf('Strain value: %.4f',isoVal),'FontSize',24);                                
            end
            
            function slideAlpha(hObj,callbackdata)
                alphaVal = get(hObj,'Value');                 
                drawIsosurface(input,alphaVal,isoVal,strain_mask);
            end 
            
            function slideIsosurfaceReal(hObj,callbackdata)
                isoVal = get(hObj,'Value');                 
                drawIsosurface(input,alphaVal,isoVal,strain_mask);
            end  
        end    
        
        function overlayWaveVectors(postPhasor, ax)
            axes(ax);
            hold on;
            scalingFactor = 0.2e-10*abs(postPhasor.plotting.object.vector1(1));
            arrowWidth = 3;
            labelsShiftRatio = 0.02;
            
            q_vector = postPhasor.experiment.qVectorLab.*scalingFactor;
            k_i = postPhasor.experiment.k_i.*scalingFactor;
            k_s = postPhasor.experiment.k_s.*scalingFactor;
            
            ki_arrow = quiver3(0,0,0,k_i(1),k_i(2),k_i(3));hold on;
            set(ki_arrow,  'Linewidth', arrowWidth, 'AutoScale', 'off','Color', 'black');
            text(k_i(1)*(1+labelsShiftRatio), k_i(2)*(1+labelsShiftRatio), k_i(3)*(1+labelsShiftRatio), 'k_{i}', 'Color', 'black', 'FontSize', 14);
            
            ks_arrow = quiver3(0,0,0,k_s(1),k_s(2),k_s(3));hold on;
            set(ks_arrow,  'Linewidth', arrowWidth, 'AutoScale', 'off','Color', 'black');
            text(k_s(1)*(1+labelsShiftRatio), k_s(2)*(1+labelsShiftRatio), k_s(3)*(1+labelsShiftRatio), 'k_{s}', 'Color', 'black', 'FontSize', 14);
                        
            q_arrow = quiver3(0,0,0,q_vector(1),q_vector(2),q_vector(3),'LineWidth',2);hold off;axis image;
            set(q_arrow,  'Linewidth', arrowWidth, 'AutoScale', 'off', 'Color', 'black');
            text(q_vector(1)*(1+labelsShiftRatio), q_vector(2)*(1+labelsShiftRatio), q_vector(3)*(1+labelsShiftRatio), 'Q', 'Color', 'black', 'FontSize', 14);            
        end
        
        function iso3dObject(postPhasor,vectorsFlag,cmap)
            % Use an input parameter to show other complex valued matrix
            if nargin == 1
                cmap = 'jet';    
                vectorsFlag = 0;
            elseif nargin == 2
                cmap = 'jet';                  
            end     
            
            input = postPhasor.object./max(abs(postPhasor.object(:)));       
            handle = figure;            
            panel = uipanel('Parent',handle,'Title','Amp_Phase','FontSize',...
            12,'Units','Normalized','Position',[.1 0.1 .77 .85]);         
            ax = axes('Parent',panel);
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.1,'Units','Normalized',...
            'Position', [0.35 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReal); 
            isoVal = 0.1;
            drawIsosurface(input,isoVal,cmap);
            
            function drawIsosurface(input,isoVal,cmap)
                cla(ax);
                axes(ax);
                isosurface(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,postPhasor.plotting.object.vector3,abs(input),isoVal,angle(input));
                
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
                
                xlim([postPhasor.plotting.object.vector1(1), postPhasor.plotting.object.vector1(end)]);
                ylim([postPhasor.plotting.object.vector2(1), postPhasor.plotting.object.vector2(end)]);
                zlim([postPhasor.plotting.object.vector3(1), postPhasor.plotting.object.vector3(end)]);
                
                rotate3d on;                
                daspect([1,1,1]);
                axis equal;
                axis vis3d;
                grid on;    
                
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];           
                colormap(cmap);
                
                if vectorsFlag
                    overlayWaveVectors(postPhasor, ax)
                end
            end
            
            function slideIsosurfaceReal(hObj,callbackdata)
                isoVal = get(hObj,'Value');                 
                drawIsosurface(input,isoVal,cmap);
            end  
        end
        
        function iso3dStrainLab(postPhasor,vectorsFlag,cmap)
            % Use an input parameter to show other complex valued matrix
            if nargin == 1
                cmap = 'jet';    
                vectorsFlag = 0;
            elseif nargin == 2
                cmap = 'jet';                  
            end                  
            
            handle = figure;            
            panel = uipanel('Parent',handle,'Title','Strain lab','FontSize',...
            12,'Units','Normalized','Position',[.1 0.1 .77 .85]);         
            ax = axes('Parent',panel);
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.1,'Units','Normalized',...
            'Position', [0.35 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReal); 
            isoVal = 0.1;
            drawIsosurface(isoVal,cmap);
            
            function drawIsosurface(isoVal,cmap)
                cla(ax);
                axes(ax);
                isosurface(postPhasor.plotting.object.vector2,...
                    postPhasor.plotting.object.vector1,...
                    postPhasor.plotting.object.vector3,...
                    abs(postPhasor.object),isoVal,postPhasor.strain);
                
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
                
                xlim([postPhasor.plotting.object.vector1(1), postPhasor.plotting.object.vector1(end)]);
                ylim([postPhasor.plotting.object.vector2(1), postPhasor.plotting.object.vector2(end)]);
                zlim([postPhasor.plotting.object.vector3(1), postPhasor.plotting.object.vector3(end)]);
                
                rotate3d on;                
                daspect([1,1,1]);
                axis equal;
                axis vis3d;
                grid on;    
%                 colorbar; 
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];           
                colormap(cmap);              
                
                if vectorsFlag
                    overlayWaveVectors(postPhasor,ax);
                end
            end
            
            function slideIsosurfaceReal(hObj,callbackdata)
                isoVal = get(hObj,'Value');  
                drawIsosurface(isoVal,cmap);
            end  
        end
        
        function iso3d_data(postPhasor,input,cmap)
            % Use an input parameter to show fft mod of the object
            if nargin == 1
                cmap = 'jet';
                input = log10(postPhasor.data)./max(log10(postPhasor.data(:)));
                title_val = 'Data';
            elseif nargin == 2
                cmap = 'jet';
                input = log10(postPhasor.object_fft_mod)./max(log10(postPhasor.object_fft_mod(:)));
                title_val = 'Reconstruction';
            end
            
            handle = figure;
            panel = uipanel('Parent',handle,'Title',title_val,'FontSize',...
            12,'Units','Normalized','Position',[.1 0.1 .77 .85]); 
        
            ax = axes('Parent',panel);
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.5,'Units','Normalized',...
            'Position', [0.1 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReal); 
            isoVal = 0.5;
            drawIsosurface(input,isoVal);
            
            function drawIsosurface(input,isoVal)
                cla(ax);
                axes(ax);
                isosurface(input,isoVal);                
                rotate3d on;
                grid on;
                axis tight;
                axis equal; 
%                 axis vis3d;
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];            
                colormap(cmap);                 
            end
            
            function slideIsosurfaceReal(hObj,callbackdata)
                isoVal = get(hObj,'Value');                 
                drawIsosurface(input,isoVal);
            end  
        end
        
        function setZoom(postPhasor,zoom_value)
            if numel(zoom_value) == 1
                postPhasor.plotting.zoom_value = repmat(zoom_value,[1,3]);
            elseif numel(zoom_value) ~= 3
                error('- Wrong number of zoom values!');
            else
                postPhasor.plotting.zoom_value = zoom_value;
                fprintf('+ Zoom values are set to: [%.2f, %.2f, %.2f]!\n',postPhasor.plotting.zoom_value(1),postPhasor.plotting.zoom_value(2),postPhasor.plotting.zoom_value(3));                
            end
        end
        
        function plot_alpha_slice(postPhasor,zoom_value,type)
            % Plot derivatives with opacity map from amplitude of the
            % reconstruction
            if nargin == 1
                zoom_value = [1,1,1];
                type = 'amp-disp';
            elseif nargin == 2
                type = 'amp-disp';
            end
            
            postPhasor.setZoom(zoom_value); % should be 3-element vector
            
            switch type
                case 'amp-disp'
                    plotLabel = 'amplitude+displacement ';
                    
                    dataImg = postPhasor.mask(:,:,round(end/2)).*postPhasor.displacement(:,:,round(end/2));
                    alphaM = abs(postPhasor.mask(:,:,round(end/2)).*postPhasor.object(:,:,round(end/2)))./max(max(abs(postPhasor.mask(:,:,round(end/2)).*postPhasor.object(:,:,round(end/2)))));
                    
                    figure('Position',[100 100 2000 500]);
                    subplot(1,3,1);imagesc(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,dataImg,'AlphaData',alphaM);                    
                    axis image;title(sprintf('%s [1,2]',plotLabel));xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));
                    colorbar;
                    
                    dataImg = squeeze(postPhasor.mask(:,round(end/2),:).*postPhasor.displacement(:,round(end/2),:));
                    alphaM = squeeze(abs(postPhasor.mask(:,round(end/2),:).*postPhasor.object(:,round(end/2),:))./max(max(abs(postPhasor.mask(:,round(end/2),:).*postPhasor.object(:,round(end/2),:)))));
                    
                    subplot(1,3,2);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector1,dataImg,'AlphaData',alphaM);
                    axis image;title(sprintf('%s [1,3]',plotLabel));xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));
                    colorbar;
                    
                    dataImg = squeeze(postPhasor.mask(round(end/2),:,:).*postPhasor.displacement(round(end/2),:,:));
                    alphaM = squeeze(abs(postPhasor.mask(round(end/2),:,:).*postPhasor.object(round(end/2),:,:))./max(max(abs(postPhasor.mask(round(end/2),:,:).*postPhasor.object(round(end/2),:,:)))));
                    
                    subplot(1,3,3);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector2,dataImg,'AlphaData',alphaM);
                    axis image;title(sprintf('%s [2,3]',plotLabel));xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
                    colorbar;
                    
                    colormap jet     
                case 'amp-strain'
                    plotLabel = 'amplitude+strain ';
                    
                    dataImg = postPhasor.mask(:,:,round(end/2)).*postPhasor.strain(:,:,round(end/2));
                    alphaM = abs(postPhasor.mask(:,:,round(end/2)).*postPhasor.object(:,:,round(end/2)))./max(max(abs(postPhasor.mask(:,:,round(end/2)).*postPhasor.object(:,:,round(end/2)))));
                    
                    figure('Position',[100 100 2000 500]);
                    subplot(1,3,1);imagesc(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,dataImg,'AlphaData',alphaM);                    
                    axis image;title(sprintf('%s [1,2]',plotLabel));xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));
                    colorbar;
                    
                    dataImg = squeeze(postPhasor.mask(:,round(end/2),:).*postPhasor.strain(:,round(end/2),:));
                    alphaM = squeeze(abs(postPhasor.mask(:,round(end/2),:).*postPhasor.object(:,round(end/2),:))./max(max(abs(postPhasor.mask(:,round(end/2),:).*postPhasor.object(:,round(end/2),:)))));
                    
                    subplot(1,3,2);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector1,dataImg,'AlphaData',alphaM);
                    axis image;title(sprintf('%s [1,3]',plotLabel));xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));
                    colorbar;
                    
                    dataImg = squeeze(postPhasor.mask(round(end/2),:,:).*postPhasor.strain(round(end/2),:,:));
                    alphaM = squeeze(abs(postPhasor.mask(round(end/2),:,:).*postPhasor.object(round(end/2),:,:))./max(max(abs(postPhasor.mask(round(end/2),:,:).*postPhasor.object(round(end/2),:,:)))));
                    
                    subplot(1,3,3);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector2,dataImg,'AlphaData',alphaM);
                    axis image;title(sprintf('%s [2,3]',plotLabel));xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
                    colorbar;
                    
                    colormap jet     
            end            
        end
        
        function plot_amplitude_slice(postPhasor,zoom_value)            
            
            postPhasor.setZoom(zoom_value); % should be 3-element vector
            
            figure('Position',[100 100 2000 500]);
            subplot(1,3,1);imagesc(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,abs(postPhasor.object(:,:,round(end/2))));
            axis image;title('Amplitude, dimensions [1,2]');xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));

            subplot(1,3,2);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector1,squeeze(abs(postPhasor.object(:,round(end/2),:))));
            axis image;title('Amplitude, dimensions [1,3]');xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));

            subplot(1,3,3);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector2,squeeze(abs(postPhasor.object(round(end/2),:,:))));
            axis image;title('Amplitude, dimensions [2,3]');xlabel('Position [nm]');ylabel('Position [nm]');
            colormap bone
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
        end
        
        function plot_phase_slice(postPhasor,zoom_value)  
            
            postPhasor.setZoom(zoom_value); % should be 3-element vector
            
            % Phase central slices
            figure('Position',[100 100 2000 500]);
            subplot(1,3,1);imagesc(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,...
                angle(postPhasor.object(:,:,round(end/2))).*postPhasor.mask(:,:,round(end/2)),'AlphaData',postPhasor.mask(:,:,round(end/2)));
            axis image;title('Phase');colorbar;xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));

            subplot(1,3,2);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector1,...
                squeeze(angle(postPhasor.object(:,round(end/2),:))).*squeeze(postPhasor.mask(:,round(end/2),:)),'AlphaData',squeeze(postPhasor.mask(:,round(end/2),:)));
            axis image;title('Phase');colorbar;xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));

            subplot(1,3,3);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector2,...
                squeeze(angle(postPhasor.object(round(end/2),:,:))).*squeeze(postPhasor.mask(round(end/2),:,:)),'AlphaData',squeeze(postPhasor.mask(round(end/2),:,:)));
            axis image;title('Phase');colorbar;xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
            colormap jet
        end
        
        function plot_displacement_slice(postPhasor, zoom_value) 
            
            postPhasor.setZoom(zoom_value); % should be 3-element vector
            
            % Phase central slices
            figure('Position',[100 100 2000 500]);
            subplot(1,3,1);imagesc(postPhasor.plotting.object.vector2,postPhasor.plotting.object.vector1,(postPhasor.displacement(:,:,round(end/2))).*postPhasor.mask(:,:,round(end/2)),'AlphaData',postPhasor.mask(:,:,round(end/2)));
            axis image;title('Displacement');colorbar;xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));
            
            subplot(1,3,2);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector1,squeeze((postPhasor.displacement(:,round(end/2),:))).*squeeze(postPhasor.mask(:,round(end/2),:)),'AlphaData',squeeze(postPhasor.mask(:,round(end/2),:)));
            axis image;title('Displacement');colorbar;xlabel('Position [nm]');ylabel('Position [nm]');
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));
            
            subplot(1,3,3);imagesc(postPhasor.plotting.object.vector3,postPhasor.plotting.object.vector2,squeeze((postPhasor.displacement(round(end/2),:,:))).*squeeze(postPhasor.mask(round(end/2),:,:)),'AlphaData',squeeze(postPhasor.mask(round(end/2),:,:)));
            axis image;title('Displacement');colorbar;xlabel('Position [nm]');ylabel('Position [nm]');
            colormap jet
            set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
        end
        
        function plot_strain_slice(postPhasor,zoom_value,slice)
            try
                fprintf('Zoom values for the plots: [%.1f,%.1f,%.1f]',zoom_value(1),zoom_value(2),zoom_value(3));
            catch
                zoom_value = [1,1,1];   
                fprintf('Zoom values are set to default: [%.1f,%.1f,%.1f]',zoom_value(1),zoom_value(2),zoom_value(3));
            end
            
            postPhasor.setZoom(zoom_value); % should be 3-element vector
            
            if nargin == 2
                try
                    figure('Position',[100 100 2000 500]);
                    subplot(1,3,1);
                    imagesc(postPhasor.plotting.strain.vector2,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            postPhasor.strain(:,:,round(end/2)).*postPhasor.strain_mask(:,:,round(end/2)),...
                            'AlphaData',postPhasor.strain_mask(:,:,round(end/2)));                    
                    axis image;colorbar;colorbar;title('Strain');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));
                    try
                        caxis(postPhasor.plotting.strain_limits(1,:))
                    catch
                    end
                        
                    subplot(1,3,2);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            squeeze(postPhasor.strain(:,round(end/2),:)).*squeeze(postPhasor.strain_mask(:,round(end/2),:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(:,round(end/2),:)));
                    axis image;colorbar;colorbar;title('Strain');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));
                    try
                        caxis(postPhasor.plotting.strain_limits(2,:))
                    catch
                    end
                    
                    subplot(1,3,3);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector2,...
                            squeeze(postPhasor.strain(round(end/2),:,:)).*squeeze(postPhasor.strain_mask(round(end/2),:,:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(round(end/2),:,:)));
                    axis image;colorbar;colorbar;title('Strain');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
                    try
                        caxis(postPhasor.plotting.strain_limits(3,:))
                    catch
                    end
                    
                    colormap jet
                catch
                    error('No strain found!')
                end
            else
                try
                    figure('Position',[100 100 2000 500]);
                    subplot(1,3,1);
                    imagesc(postPhasor.plotting.strain.vector2,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            postPhasor.strain(:,:,slice(3)).*postPhasor.strain_mask(:,:,slice(3)),...
                            'AlphaData',postPhasor.strain_mask(:,:,slice(3)));                    
                    axis image;colorbar;colorbar;title('Strain');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));

                    subplot(1,3,2);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            squeeze(postPhasor.strain(:,slice(2),:)).*squeeze(postPhasor.strain_mask(:,slice(2),:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(:,slice(2),:)));
                    axis image;colorbar;colorbar;title('Strain');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));

                    subplot(1,3,3);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector2,...
                            squeeze(postPhasor.strain(slice(1),:,:)).*squeeze(postPhasor.strain_mask(slice(1),:,:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(slice(1),:,:)));
                    axis image;colorbar;colorbar;title('Strain');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
                    colormap jet
                catch
                    error('No strain found!')
                end
            end
        end
        
        function plotStrainSegmentsSlice(postPhasor,zoom_value)
           postPhasor.setZoom(zoom_value); % should be 3-element vector        
                try
                    % Core
                    figure('Position',[100 100 2000 500]);
                    subplot(1,3,1);
                    imagesc(postPhasor.plotting.strain.vector2,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            postPhasor.strain(:,:,round(end/2)).*postPhasor.strain_mask_bulk(:,:,round(end/2)),...
                            'AlphaData',postPhasor.strain_mask(:,:,round(end/2)));                    
                    axis image;colorbar;colorbar;title('Strain: bulk');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));

                    subplot(1,3,2);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            squeeze(postPhasor.strain(:,round(end/2),:)).*squeeze(postPhasor.strain_mask_bulk(:,round(end/2),:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(:,round(end/2),:)));
                    axis image;colorbar;colorbar;title('Strain: bulk');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));

                    subplot(1,3,3);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector2,...
                            squeeze(postPhasor.strain(round(end/2),:,:)).*squeeze(postPhasor.strain_mask_bulk(round(end/2),:,:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(round(end/2),:,:)));
                    axis image;colorbar;colorbar;title('Strain: bulk');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
                    colormap jet
                    
                    % Shell
                    figure('Position',[100 100 2000 500]);
                    subplot(1,3,1);
                    imagesc(postPhasor.plotting.strain.vector2,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            postPhasor.strain(:,:,round(end/2)).*postPhasor.strain_mask_shell(:,:,round(end/2)),...
                            'AlphaData',postPhasor.strain_mask(:,:,round(end/2)));                    
                    axis image;colorbar;colorbar;title('Strain: shell');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(1));

                    subplot(1,3,2);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector1(2:end),...
                            squeeze(postPhasor.strain(:,round(end/2),:)).*squeeze(postPhasor.strain_mask_shell(:,round(end/2),:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(:,round(end/2),:)));
                    axis image;colorbar;colorbar;title('Strain: shell');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(2));

                    subplot(1,3,3);
                    imagesc(postPhasor.plotting.strain.vector3,...
                            postPhasor.plotting.strain.vector2,...
                            squeeze(postPhasor.strain(round(end/2),:,:)).*squeeze(postPhasor.strain_mask_shell(round(end/2),:,:)),...
                            'AlphaData',squeeze(postPhasor.strain_mask(round(end/2),:,:)));
                    axis image;colorbar;colorbar;title('Strain: shell');xlabel('Position [nm]');ylabel('Position [nm]');
                    set(gca,'FontSize',20);zoom(postPhasor.plotting.zoom_value(3));
                    colormap jet
                catch
                    error('No strain or core-shell mask found!')
                end    
        end        
        
        function plot_central_profiles(postPhasor,volume)
            if ~exist('volume','var')
                volume = 'all';
            end
            
            if ismember(volume,['amplitude', 'all'])
                figure('Position',[100 100 2000 500]);
                subplot(1,3,1);
                val = squeeze(abs(postPhasor.object(:,round(end/2),round(end/2))));
                val_contour = squeeze((postPhasor.mask(:,round(end/2),round(end/2)))).*max(val(:));
                val_contour_mask = val_contour>0;
                val_contour_length = sum(val_contour_mask).*postPhasor.object_sampling(1);
                
                plot(postPhasor.plotting.object.vector1,val,'.-');hold on
                plot(postPhasor.plotting.object.vector1,val_contour,'--');
                title(sprintf('Amplitude profile: %.2f nm',val_contour_length*1e9));
                grid on;xlabel('Position y [nm]');set(gca,'FontSize',20);

                subplot(1,3,2);
                val = squeeze(abs(postPhasor.object(round(end/2),:,round(end/2))));
                val_contour = squeeze((postPhasor.mask(round(end/2),:,round(end/2)))).*max(val(:));                
                val_contour_mask = val_contour>0;
                val_contour_length = sum(val_contour_mask).*postPhasor.object_sampling(2);
                plot(postPhasor.plotting.object.vector2,val,'.-');hold on
                plot(postPhasor.plotting.object.vector2,val_contour,'--');
                title(sprintf('Amplitude profile: %.2f nm',val_contour_length*1e9));
                grid on;xlabel('Position x [nm]');set(gca,'FontSize',20);

                subplot(1,3,3);
                val = squeeze(abs(postPhasor.object(round(end/2),round(end/2),:)));
                val_contour = squeeze(postPhasor.mask(round(end/2),round(end/2),:)).*max(val(:));
                val_contour_mask = val_contour>0;
                val_contour_length = sum(val_contour_mask).*postPhasor.object_sampling(3);
                plot(postPhasor.plotting.object.vector3,val,'.-');hold on
                plot(postPhasor.plotting.object.vector3,val_contour,'--');
                title(sprintf('Amplitude profile: %.2f nm',val_contour_length*1e9));
                grid on;xlabel('Position z [nm]');set(gca,'FontSize',20);
            end
            if ismember(volume,['phase', 'all'])       
                figure;
                subplot(1,3,1);
                val = squeeze(angle(postPhasor.object(:,round(end/2),round(end/2))));
                plot(postPhasor.plotting.object.vector1,val,'.-');hold on
                plot(postPhasor.plotting.object.vector1,squeeze((postPhasor.mask(:,round(end/2),round(end/2)))).*max(val(:)),'--');
                title('Phase profile');grid on;xlabel('Position [nm]');

                subplot(1,3,2);
                val = squeeze(angle(postPhasor.object(round(end/2),:,round(end/2))));
                plot(postPhasor.plotting.object.vector2,val,'.-');hold on
                plot(postPhasor.plotting.object.vector2,squeeze((postPhasor.mask(round(end/2),:,round(end/2)))).*max(val(:)),'--');
                title('Phase profile');grid on;xlabel('Position [nm]');

                subplot(1,3,3);
                val = squeeze(angle(postPhasor.object(round(end/2),round(end/2),:)));
                plot(postPhasor.plotting.object.vector3,val,'.-');hold on
                plot(postPhasor.plotting.object.vector3,squeeze(postPhasor.mask(round(end/2),round(end/2),:)).*max(val(:)),'--');
                title('Phase profile');grid on;xlabel('Position [nm]');
            end
            if ismember(volume,['displacement', 'all'])            
                figure;
                subplot(1,3,1);
                val = squeeze(postPhasor.displacement(:,round(end/2),round(end/2)));
                plot(postPhasor.plotting.object.vector1,val,'.-');hold on
                plot(postPhasor.plotting.object.vector1,squeeze((postPhasor.mask(:,round(end/2),round(end/2)))).*max(val(:)),'--');
                title('Displacement profile');grid on;xlabel('Position [nm]');

                subplot(1,3,2);
                val = squeeze(postPhasor.displacement(round(end/2),:,round(end/2)));
                plot(postPhasor.plotting.object.vector2,val,'.-');hold on
                plot(postPhasor.plotting.object.vector2,squeeze((postPhasor.mask(round(end/2),:,round(end/2)))).*max(val(:)),'--');
                title('Displacement profile');grid on;xlabel('Position [nm]');

                subplot(1,3,3);
                val = squeeze(postPhasor.displacement(round(end/2),round(end/2),:));
                plot(postPhasor.plotting.object.vector3,val,'.-');hold on
                plot(postPhasor.plotting.object.vector3,squeeze(postPhasor.mask(round(end/2),round(end/2),:)).*max(val(:)),'--');
                title('Displacement profile');grid on;xlabel('Position [nm]');
            end
            if ismember(volume,['strain', 'all'])      
                figure;
                subplot(1,3,1);
                val = squeeze(postPhasor.strain(:,round(end/2),round(end/2)));
                plot(postPhasor.plotting.strain.vector1,val,'.-');hold on
                plot(postPhasor.plotting.strain.vector1,squeeze((postPhasor.strain_mask(:,round(end/2),round(end/2)))).*max(val(:)),'--');
                title('Strain profile');grid on;xlabel('Position [nm]');

                subplot(1,3,2);
                val = squeeze(postPhasor.strain(round(end/2),:,round(end/2)));
                plot(postPhasor.plotting.strain.vector2,val,'.-');hold on
                plot(postPhasor.plotting.strain.vector2,squeeze((postPhasor.strain_mask(round(end/2),:,round(end/2)))).*max(val(:)),'--');
                title('Strain profile');grid on;xlabel('Position [nm]');

                subplot(1,3,3);
                val = squeeze(postPhasor.strain(round(end/2),round(end/2),:));
                plot(postPhasor.plotting.strain.vector3,val,'.-');hold on
                plot(postPhasor.plotting.strain.vector3,squeeze(postPhasor.strain_mask(round(end/2),round(end/2),:)).*max(val(:)),'--');
                title('Strain profile');grid on;xlabel('Position [nm]');    
            end            
        end
        
        % Save functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save the current instance of the object
        function save_instance(postPhasor)
            postPhasor.dataTime = getTimeStamp;  
            % Save the whole instance
            save_path = fullfile(postPhasor.pre_path, 'post_processing');
            mkdir(save_path);
            s = fullfile(save_path, 'postPhasor.mat');
            save(s,'postPhasor');
            fprintf('+ postPhasor instance is saved to %s\n', s);
        end
        
        function save_figures(postPhasor, zoom_value, values)            
            try 
                zoom_value;                
            catch
                zoom_value = [1,1,1];
                warning('No zoom value was specified! Using [1,1,1]...');
            end
            
            try
                values;
            catch
                values = {'all'};
                warning('No values were specified for plotting! Using all...');
            end
            
            save_path = fullfile(postPhasor.pre_path, 'post_processing');
            save_path_figures = fullfile(save_path, 'figures');
            
            mkdir(save_path);
            mkdir(save_path_figures);
            
            if ismember('amp',values) || ismember('all',values)
                % Figures
                postPhasor.plot_amplitude_slice(zoom_value);
                hFig = gcf;
                set(hFig,'Units','Inches');
                pos = get(hFig,'Position');
                set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
                print(hFig,fullfile(save_path_figures,'amplitude_slice.pdf'),'-dpdf','-r0');
                print(hFig,fullfile(save_path_figures,'amplitude_slice.png'),'-dpng','-r0');
                savefig(hFig,fullfile(save_path_figures,'amplitude_slice.fig'));
                close(hFig);
            
            elseif ismember('phase',values) || ismember('all',values)
                postPhasor.plot_phase_slice(zoom_value);
                hFig = gcf;
                set(hFig,'Units','Inches');
                pos = get(hFig,'Position');
                set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
                print(hFig,fullfile(save_path_figures,'phase_slice.pdf'),'-dpdf','-r0');
                print(hFig,fullfile(save_path_figures,'phase_slice.png'),'-dpng','-r0');
                savefig(hFig,fullfile(save_path_figures,'phase_slice.fig'));
                close(hFig);
            
            elseif ismember('strain',values) || ismember('all',values)
                postPhasor.plot_strain_slice(zoom_value);
                hFig = gcf;
                set(hFig,'Units','Inches');
                pos = get(hFig,'Position');
                set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
                print(hFig,fullfile(save_path_figures,'strain_slice.pdf'),'-dpdf','-r0');
                print(hFig,fullfile(save_path_figures,'strain_slice.png'),'-dpng','-r0');
                savefig(hFig,fullfile(save_path_figures,'strain_slice.fig'));
                close(hFig);
            
            elseif ismember('displ',values) || ismember('all',values)
                postPhasor.plot_displacement_slice(zoom_value);
                hFig = gcf;
                set(hFig,'Units','Inches');
                pos = get(hFig,'Position');
                set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
                print(hFig,fullfile(save_path_figures,'displacement_slice.pdf'),'-dpdf','-r0');
                print(hFig,fullfile(save_path_figures,'displacement_slice.png'),'-dpng','-r0');
                savefig(hFig,fullfile(save_path_figures,'displacement_slice.fig'));
                close(hFig);
            
            elseif ismember('amp_profiles',values) || ismember('all',values)
                postPhasor.plot_central_profiles('amplitude');
                hFig = gcf;
                set(hFig,'Units','Inches');
                pos = get(hFig,'Position');
                set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
                print(hFig,fullfile(save_path_figures,'amplitude_central_profiles.pdf'),'-dpdf','-r0');
                print(hFig,fullfile(save_path_figures,'amplitude_central_profiles.png'),'-dpng','-r0');
                savefig(hFig,fullfile(save_path_figures,'amplitude_central_profiles.fig'));
                close(hFig);
            end
            
            fprintf('+ Figures are saved to %s\n', save_path_figures);
        end
        
        function save_vtk(postPhasor)
             save_path = fullfile(postPhasor.pre_path, 'post_processing');
            mkdir(save_path);
            
            savemat2vtk(fullfile(save_path, 'object.vtk'),          postPhasor.object.*postPhasor.mask,     postPhasor.object_sampling*1e9);
            try
                savemat2vtk(fullfile(save_path, 'displacement.vtk'), abs(postPhasor.object), postPhasor.object_sampling*1e9,     postPhasor.displacement.*postPhasor.mask);
            catch
            end
            try
                savemat2vtk(fullfile(save_path, 'strain.vtk'),      abs(postPhasor.object), postPhasor.object_sampling*1e9,     postPhasor.strain.*postPhasor.strain_mask);
            catch
            end
            
            try
                scalingFactor = 0.2e-10*abs(postPhasor.plotting.object.vector1(1));                        
                q_vector = postPhasor.experiment.qVectorLab.*scalingFactor;
                saveVec2Vtk(fullfile(save_path, 'qVectorLab.vtk'),[0,0,0],q_vector);
            catch
            end
            
            fprintf('+ VTK files are saved to %s\n', save_path);
        end
        
        function save_all(postPhasor,zoom_value)      
            postPhasor.dataTime = getTimeStamp;                        
            save_path = fullfile(postPhasor.pre_path, 'post_processing');
            save_path_figures = fullfile(save_path, 'figures');
            
            mkdir(save_path);
            mkdir(save_path_figures);
            
            postPhasor.save_instance;
            postPhasor.save_figures(zoom_value); 
            postPhasor.save_vtk;
            
            postPhasor.segment_strain_mask;
            postPhasor.calculate_strain_histogram('bulk');
            
            hFig = gcf;
            set(hFig,'Units','Inches');
            pos = get(hFig,'Position');
            set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
            print(hFig,fullfile(save_path_figures,'histogram_bulk.pdf'),'-dpdf','-r0');
            print(hFig,fullfile(save_path_figures,'histogram_bulk.png'),'-dpng','-r0');
            savefig(hFig,fullfile(save_path_figures,'histogram_bulk.fig'));
            close(hFig);
            
            postPhasor.calculate_strain_histogram('shell');
            hFig = gcf;
            set(hFig,'Units','Inches');
            pos = get(hFig,'Position');
            set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
            print(hFig,fullfile(save_path_figures,'histogram_shell.pdf'),'-dpdf','-r0');
            print(hFig,fullfile(save_path_figures,'histogram_shell.png'),'-dpng','-r0');
            savefig(hFig,fullfile(save_path_figures,'histogram_shell.fig'));
            close(hFig);
            
            fprintf('Saved everything to: %s\n',save_path);
        end
                
    end
end