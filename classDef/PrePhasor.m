classdef PrePhasor < handle
% Naming conventions

% % Methods
% Method names should be all lower case
% Words in an method name should be separated by an underscore
% Non-public method should begin with a single underscore
% If a method name needs to be mangled, two underscores may begin its name

% % Constants
% Constant names must be fully capitalized
% Words in a constant name should be separated by an underscore
 
% % Instance Variables
% Instance variable names should be all lower case
% Words in an instance variable name should be separated by an underscore
% Non-public instance variables should begin with a single underscore
% If an instance name needs to be mangled, two underscores may begin its name
    
% Structure of data field: [vertical_axis | horizontal_axis | angle/energy_axis | scan_axis]

    properties
        data;
        mask;    
        crop_flag = 0;
        white_field;
        dark_field;
        data_max;
        data_min;
        data_average;
        data_binned;
        data_cropped;
        data_3d;
        data_integral;
        data_meta;           
    end
    
    % - Constructor of the object -
    methods
        function prePhasor = PrePhasor(input_param)          
            prePhasor.data_meta      = input_param;                              
        end
    end    
    
    % - Processing methods -
    methods
        function create_path(prePhasor,path)
            switch prePhasor.data_meta.beamline_id
                
                case '34idc'
                    % For 34 idc output is a list of tiff files
                    % =========================================================================
                    % --- Assemble 'master' file names
                    % =========================================================================
                    if ~exist('path')                
                        try % 34IDC option        
                           main_dir =  fullfile(prePhasor.data_meta.pre_path_data,prePhasor.data_meta.beamtime_id,...
                                                ['AD' prePhasor.data_meta.beamtime_prefix '_' prePhasor.data_meta.sample_name]);
                        catch
                           main_dir = fullfile('/asap3','petra3','gpfs','p10','2018','data', prePhasor.data_meta.beamtime_id, 'raw'); % E.g. on Office Linux PCs            ; % E.g. on  Windows PCs        
                        end
                    else
                        main_dir = path;
                    end
                    
                    if strcmpi(prePhasor.data_meta.sample_name(end),'_') ~= 1
                        prePhasor.data_meta.sample_name = [prePhasor.data_meta.sample_name '_'];
                    end

                    for jj = 1:numel(prePhasor.data_meta.scan_number)
                        master_folder = fullfile(main_dir, [prePhasor.data_meta.beamtime_prefix ...
                                                 prePhasor.data_meta.sample_name 'S' sprintf('%04i',prePhasor.data_meta.scan_number(jj))]);

                        t = dir(fullfile(master_folder,'*.tif'));

                        for ii = 1:length(t)
                            prePhasor.data_meta.scan(jj).file(ii).name        = fullfile(master_folder,t(ii).name);
                        end
                    end
                    
                case 'nanomax'
                    if ~exist('path')                	         
                        try % NanoMAx option                            
                            main_dir =  fullfile(prePhasor.data_meta.pre_path_data,'MAXIV',prePhasor.data_meta.beamline_id,prePhasor.data_meta.beamtime_id,'raw',...
                                                 prePhasor.data_meta.sample_name);
                        catch
                            warning('no path!');
                        end
                    else
                        main_dir = path;
                    end
                    
                    for jj = 1:numel(prePhasor.data_meta.scan_number)                        
                        prePhasor.data_meta.scan(jj).file.name        = fullfile(main_dir,['scan_' sprintf('%04i_',prePhasor.data_meta.scan_number(jj)) prePhasor.data_meta.detector_id '_0000.hdf5']);                        
                    end
            end
        end
        
        function create_mask(prePhasor,dim)
            function f_capturekeystroke(H,E)
                disp(E.Key);
                switch E.Key
                    case 'escape'
                        fprintf('Mask creation is broken at:\n %s\n',[prePhasor.data_meta.sample_name ' | Scan '...
                                       num2str(prePhasor.data_meta.scan_number(jj)) ' | Frame ' ...
                                       num2str(ii)]);                                                         
                        flag_exit       = 1;  
                    case 'space'
                        flag_next_frame = 1;
                        disp('Frame skipped!');
                    case 'control'
                        flag_control = 1;
                        disp('Frame masking!');
                end
            end
            
            if nargin==1
                dim = '3D';
            end
            switch dim
                case '2D'
                    for jj = 1:size(prePhasor.data,4)
                        jj
                    end
                case '3D'   
                    flag_exit       = 0;
                    hF = figure('keypressfcn',@f_capturekeystroke);
                    hAx = axes('Parent',hF);
                    disp('Masking:\n esc - abort;\n space - next frame;\n Ctrl - mask frame;\n')
                    for jj = 1:size(prePhasor.data,4)
                        if flag_exit
                            return;
                        else
                            for ii = 1:size(prePhasor.data,3)   
                                if flag_exit
                                    return;
                                else
                                    cla(hAx);
                                    flag_next_frame = 0;  
                                    flag_control    = 0;
                                    prePhasor.mask(:,:,ii,jj) = zeros(size(prePhasor.data(:,:,ii,jj)));
                                    while ~flag_next_frame & ~flag_exit                                                                                
                                        imagesc((prePhasor.data(:,:,ii,jj)),[0 1]);
                                        axis image;
                                        colormap hot;
                                        colormap jet;
                                        title({[prePhasor.data_meta.sample_name ' | Scan '...
                                               num2str(prePhasor.data_meta.scan_number(jj)) ' | Frame ' ...
                                               num2str(ii)], 'Space - next frame | Ctrl - mask | Esc - exit'});
                                        if flag_exit
                                            close(hF);
                                            return;                                            
                                        else
                                            waitforbuttonpress;
                                            if flag_exit
                                                close(hF);
                                                return;
                                            else
                                                if flag_control
                                                    hROI = drawfreehand(hAx);
                                                    prePhasor.mask(:,:,ii,jj) = prePhasor.mask(:,:,ii,jj)+createMask(hROI);
                                                    waitforbuttonpress; 
                                                else
                                                    disp('skipped')
                                                    break;
                                                end
                                            end
                                        end

                                                                       
                                    end 
                                    prePhasor.mask(:,:,ii,jj) = prePhasor.mask(:,:,ii,jj)>0;
                                    disp('Mask frame recorded!');
                                end
                            end
                        end
                    end
                prePhasor.mask = abs(prePhasor.mask-1);
                disp('Full 3D mask recorded!'); 
                close(hF);
            end
        end
        
        function read_tif(prePhasor)
            try
                for jj = 1:numel(prePhasor.data_meta.scan_number)
                    % Read data from 
                    for ii = 1:length(prePhasor.data_meta.scan(jj).file)
                        prePhasor.data(:,:,ii,jj) = single(imread(prePhasor.data_meta.scan(jj).file(ii).name));
                    end
                    fprintf('Loaded: %s \n',[prePhasor.data_meta.sample_name 'S' sprintf('%04i',prePhasor.data_meta.scan_number(jj))])
                end
            catch
                error('Can not load the data!')
            end
        end   
        
        function read_nanomax_merlin(prePhasor)            
            try
                % Extract scan information first                
                try                
                    for kk = 1:numel(prePhasor.data_meta.scan_number)  
                        if prePhasor.data_meta.crop_flag
                            prePhasor.data = openmultimerlin_roi(prePhasor.data_meta.scan(kk).file.name,prePhasor.data_meta.start_row,prePhasor.data_meta.end_row,...
                                [prePhasor.data_meta.roi(1),prePhasor.data_meta.roi(2),prePhasor.data_meta.roi(3),prePhasor.data_meta.roi(4),prePhasor.data_meta.start_column,prePhasor.data_meta.end_column]);
                        else
                            prePhasor.data = openmultimerlin_roi(prePhasor.data_meta.scan(kk).file.name);
                        end
                        fprintf('Loaded: %d \n',kk)
                    end
                catch
                    error('No master file!')
                end
            catch
                error('Can not load the data!')
            end
        end
        
        function read_mask(prePhasor)
%             file_temp = fullfile(scan.data_meta.save_folder,[scan.data_meta.sample_name,'_',num2str(scan.data_meta.scan_number)],scan.data_meta.mask_name);
            try
                switch prePhasor.data_meta.mask_path(end-2:end)
                    case 'mat'
                        load(prePhasor.data_meta.mask_path);
                        prePhasor.mask = mask;
                    case 'tif'
                        prePhasor.mask = single(imread(prePhasor.data_meta.mask_path));
                end
                fprintf('Mask loaded:\n %s\n',prePhasor.data_meta.mask_path);
            catch
                warning('No mask specified!');
            end
        end
               
        function read_white_field(prePhasor)
            try
                switch prePhasor.data_meta.white_field_path(end-2:end)
                    case 'mat'
                        load(prePhasor.data_meta.white_field_path);
                        prePhasor.white_field = white_field;
                    case 'tif'
                        prePhasor.white_field = single(imread(prePhasor.data_meta.white_field_path));
                end
%                 scan.white_field(scan.white_field < 6000) = 1e25;
                disp('### White field loaded ###');
            catch
                warning('No white field specified!');
            end
        end
        
        function read_dark_field(prePhasor)
            try
                switch prePhasor.data_meta.white_field_path(end-2:end)
                    case 'mat'
                        load(prePhasor.data_meta.dark_field_path);
                        prePhasor.dark_field = dark_field;
                    case 'tif'
                        prePhasor.dark_field = single(imread(prePhasor.data_meta.dark_field_path));
                end
                disp('### Dark field loaded ###');
            catch
                warning('No dark field specified!');
            end
        end        
        
        function read_beam_current(prePhasor)
             switch prePhasor.data_meta.beamline_id                
                case 'nanomax'
                    prePhasor.data_meta.nanomax.beam_current = h5read(prePhasor.data_meta.master_file_nanomax,sprintf('/entry%d/measurement/beam_current/',prePhasor.data_meta.scan_number));
             end
        end
        
        function correct_low_cutoff(prePhasor)
            prePhasor.data(prePhasor.data<=prePhasor.data_meta.low_cutoff) = 0;
            disp('Data was low-tresholded');
        end
        
        function correct_dark_field(prePhasor)
            disp('### Correcting by dark-field ###');
            try
                if ~isempty(prePhasor.dark_field) & size(prePhasor.data(:,:,1))==size(prePhasor.dark_field)
                    for jj = 1:size(prePhasor.data,4) 
                        for ii = 1:size(prePhasor.data,3) 
                            t = prePhasor.data(:,:,ii,jj);
                            t(prePhasor.dark_field>1) = 0;
                            prePhasor.data(:,:,ii,jj) = t;                                
                        end
                        fprintf('Processign Scan #%d\n',jj);
                    end
                    disp('Data corrected by dark field!')
                elseif ~isempty(prePhasor.dark_field) & size(prePhasor.data(:,:,1))~=size(prePhasor.dark_field)
                    error('Dark field size does not match data size!')
                elseif isempty(prePhasor.dark_field)
                    error('No dark field!')
                end
            catch
                if ndims(prePhasor.data) ~= 3
                    warning('The data is not 3D! Skipped. Add functionality to method!')
                else
                    error('Can not correct by dark field');
                end
            end
        end
               
        function correct_white_field(prePhasor)
            disp('### Correcting by white-field ###');
            try
                if ~isempty(prePhasor.white_field) & size(prePhasor.data(:,:,1))==size(prePhasor.white_field)                    
                    for jj = 1:size(prePhasor.data,4) 
                        for ii = 1:size(prePhasor.data,3) 
                            prePhasor.data(:,:,ii,jj) = max(prePhasor.white_field(:)).*prePhasor.data(:,:,ii,jj)./prePhasor.white_field; 
                            prePhasor.data(isinf(prePhasor.data)) = 0;
                            prePhasor.data(isnan(prePhasor.data)) = 0;
                        end
                        fprintf('Processing Scan #%d\n',jj);
                    end
                    disp('Data corrected by white field!')
                elseif ~isempty(prePhasor.white_field) & size(prePhasor.data(:,:,1))~=size(prePhasor.white_field)
                    error('White field size does not match data size!')
                elseif isempty(prePhasor.white_field)
                    error('No White field!')
                end
            catch
                if ndims(prePhasor.data) ~= 3
                    warning('The data is not 3D! Skipped. Add functionality to method!')
                else
                    error('Can not correct by white field');
                end
            end
        end
        
        function correct_mask(prePhasor)
            disp('### Masking the data ###');
            try
                prePhasor.data = prePhasor.data.*prePhasor.mask;
                disp('Mask applied!');
            catch
                warning('No mask specified or exists!');
            end
        end
        
        function correct_mask_nanomax(prePhasor)
            disp('### Masking the data ###');
            try
                for ii =1:size(prePhasor.data,3)
                    for jj =1:size(prePhasor.data,4)
                        prePhasor.data(:,:,ii,jj) = prePhasor.data(:,:,ii,jj).*single(prePhasor.mask);
                    end
                end
                disp('Mask applied!');
            catch
                warning('No mask specified or exists!');
            end
        end
        
        function crop(prePhasor)
            show_data_average(prePhasor,'log');
            hRect = drawrectangle;            
            disp('### Cropping the dataset ###')
            for jj = 1:size(prePhasor.data,4)
                for ii = 1:size(prePhasor.data,3)
                    prePhasor.data_crop(:,:,ii,jj) = imcrop(squeeze(prePhasor.data(:,:,ii,jj)),hRect.Position);
                end
            end
            figure; imagesc(log10(squeeze(mean(mean(prePhasor.data_crop,4),3)))); colormap jet; axis image;
            prePhasor.crop_flag = 1;
        end
        
        function crop_auto(prePhasor)
            % Find the center of mass
            com = ndimCOM(prePhasor.data,'auto');
            % Minimum windows in each dimension
            window = floor(min(com,size(prePhasor.data)-com));
            window(1:2) = min(window(1),window(2));
            prePhasor.data_cropped = prePhasor.data(com(1)-window(1)+1:com(1)+window(1),...
                                            com(2)-window(2)+1:com(2)+window(2),...
                                            com(3)-window(3)+1:com(3)+window(3));
        end
        
        function correct_hot_pixel(prePhasor,x,y,interpolate)
            disp('### Hot-pixels correction ###');
            if ~interpolate
                prePhasor.data(x,y,:,:) = 0;
                fprintf('Hot pixel [x:%d y:%d] zeroed!\n',x,y);
            else               
                for jj = 1:size(prePhasor.data,4)
                    for ii = 1:size(prePhasor.data,3)
                        try
                            prePhasor.data(y,x,ii,jj) = mean(mean(prePhasor.data(y-1:2:y+1,x-1:2:x+1,ii,jj)));
                        catch
                            try
                                prePhasor.data(y,x,ii,jj) = mean(mean(prePhasor.data(y-1:2:y+1,x+1,ii,jj)));                        
                            catch
                                prePhasor.data(y,x,ii,jj) = mean(mean(prePhasor.data(y+1,x-1:2:x+1,ii,jj)));                        
                            end                            
                        end
                    end
                end                                
                fprintf('Hot pixel [x:%d y:%d] interpolated!\n',x,y);
            end                
        end        
        
        function average(prePhasor,dimsAverage)
            disp('### Averaging the data ###');            
            try  
%                 clear data_average;
                prePhasor.data_average = squeeze(mean(prePhasor.data,dimsAverage));                                  
            catch
                disp('Data not-averaged!');
            end                
        end
        
        function flip(prePhasor,dim)            
            if ~exist('dim')
                warning('No dimension specified: lr / ud. Skip.')
            else
                fprintf('### Flipping the data %s ###', dim);
                switch dim
                    case 'lr'
                        prePhasor.data = fliplr(prePhasor.data);
                    case 'ud'
                        prePhasor.data = flipud(prePhasor.data);
                end
            end
        end
        
        function bin2D(prePhasor, binning_size)
            for ii = 1:size(prePhasor.data,3)
                for jj = 1:size(prePhasor.data,4)
                    convoluted = conv2(prePhasor.data(:,:,ii,jj), ones(binning_size));
                    convoluted_size = size(convoluted);
                    prePhasor.data_binned(:,:,ii,jj) = convoluted(binning_size:binning_size:convoluted_size(1), binning_size:binning_size:convoluted_size(2));
                end
            end            
        end
                
        function [COM] = ndimCOM(IN,type)    
            disp('### Calculating the center of mass ###');
            if strcmp(type,'manual')        
                imagesc(log10(IN));axis image
                h = impoly;
                hMask = createMask(h);
                IN = IN.*hMask;
            end
            C = cellfun(@(n) 1:n, num2cell(size(IN)),'uniformoutput',0);
            [C{:}] = ndgrid(C{:});
            C = cellfun(@(x) x(:), C,'uniformoutput',0);
            C = [C{:}];
            COM(:,:) = IN(:).'*C/sum(IN(:),'double');
        end        
        
        function combine(prePhasor)
            disp('### Aligning multiple-scan data with respect to 1st in array ###');
            % Center
            reference_object = 1;
            
            for jj = 1:size(prePhasor.data,4)-1
%                 try
%                     c = convn(gpuArray(scan.data(:,:,:,reference_object)),gpuArray(scan.data(end:-1:1,end:-1:1,end:-1:1,reference_object+jj)));
%                     c = gather(c);
%                 catch
                    disp('GPU acceleration was not found or too big array, therefore wait =)');                    
                    c = convn((prePhasor.data(:,:,:,reference_object)),(prePhasor.data(end:-1:1,end:-1:1,end:-1:1,reference_object+jj)));
%                 end
                
                [x,y,z] = ind2sub(size(c),find(c == max(max(max(c)))));

                shift = [x-size(prePhasor.data(:,:,:,reference_object),1),y-size(prePhasor.data(:,:,:,reference_object),2),z-size(prePhasor.data(:,:,:,reference_object),3)];

                prePhasor.data(:,:,:,reference_object+jj) = circshift(prePhasor.data(:,:,:,reference_object+jj),shift); 
                
                fprintf('Shifted Scan #%d\n', jj);
            end           
        end
        
        function normalize_exposure(prePhasor)
            disp('### Normalizing the data by exposure time ###');
            if ~isempty(prePhasor.data_meta.exposure)
                prePhasor.data = prePhasor.data./prePhasor.data_meta.exposure;
                fprintf('Data is normalized by exposure: %.3f s\n', prePhasor.data_meta.exposure);
                if isempty(prePhasor.data_meta.dead_time)
                    prePhasor.data_meta.dead_time = 0 ;
                end
                prePhasor.data = prePhasor.data./(1-prePhasor.data_meta.dead_time.*prePhasor.data);
                fprintf('Data is normalized by dead time: %.3e s\n', prePhasor.data_meta.dead_time);                
            else
                error('Exposure time is missing. Skipping...')
            end
        end
        
        function prepare_3d(prePhasor)
            try
                prePhasor.data_3d = log10(prePhasor.data)./max(max(max(log10(prePhasor.data))));
            catch
                warning('Can not prepare 3D array, checl the method!')
            end
        end
                        
        function integrate(prePhasor,dimsIntegrate) % [dimensions to integrate] 
            disp('### Integrating the data ###');            
            try
                prePhasor.data_integral = squeeze(sum(prePhasor.data,dimsIntegrate));
                fprintf('Dimensions integrated: \n %d \n',dimsIntegrate);
            catch
                disp('Data not-integrated!');
            end
        end 
        
        function maximize(prePhasor)
            disp('### Getting max values from each frame ###');
            try
                if ndims(prePhasor.data) == 4                
                    prePhasor.data_max = squeeze(sum(sum(squeeze(sum(prePhasor.data,3)),1),2));
                    disp('Data integrated sown to 3D!');
                else
                    for ii = 1:size(prePhasor.data,3)
                        prePhasor.data_max(ii) = squeeze(max(max(prePhasor.data(:,:,ii))));                        
                    end
                end                
            catch
                disp('Data not-maximised!');
            end
        end
        
        function minimize(prePhasor)
            disp('### Getting min values from each frame ###');
            try
                if ndims(prePhasor.data) == 4                
                    prePhasor.data_max = squeeze(sum(sum(squeeze(sum(prePhasor.data,3)),1),2));
                    disp('Data integrated sown to 3D!');
                else
                    for ii = 1:size(prePhasor.data,3)                        
                        m = double(prePhasor.data(:,:,ii)>0);
                        m(m==0) = NaN;
                        prePhasor.data_min(ii) = squeeze(nanmin(nanmin(m.*prePhasor.data(:,:,ii))));                        
                    end
                end                
            catch
                disp('Data not-minimized!');
            end
        end
    end
    
    % Show methods
    methods  
        function show_data_max(prePhasor)
            maximize(prePhasor);
            figure; plot(log10(prePhasor.data_max),'LineWidth',2,'Marker','o');
        end
        
        function show_data_min(prePhasor)
            minimize(prePhasor);
            figure; plot((prePhasor.data_min),'LineWidth',2,'Marker','o');
        end
        
        function show_dark_field(prePhasor)
            try
                figure;            
                imagesc(prePhasor.dark_field);
                axis image;
                colormap jet;
                colorbar;
                title('Dark field');
            catch
                error('No dark field!')
            end
        end
        
        function show_white_field(prePhasor)
            try
                figure;            
                imagesc(prePhasor.white_field);axis image;colormap jet;colorbar
                title('White field');
            catch
                error('No white field!')
            end
        end
        
        function show_3d(prePhasor,isoVal)                    
            if nargin == 1
                isoVal = 0.5;
            end                            
            try
                isosurface((prePhasor.data_3d),isoVal); axis image
            catch
                prepare_3d(prePhasor);
                isosurface((prePhasor.data_3d),isoVal); axis image
            end
        end
        
        function show_data_scroll(prePhasor,scale,max_val)
            if ~exist('scale')
                scale = 'log';
            end
            
            if nargin == 2
                average(prePhasor);
                max_val = mean(prePhasor.data_average(:))*.5;
            end
            
            if ndims(prePhasor.data) == 4    
                switch scale
                    case 'log'
                        handle = implay(log10(sum(prePhasor.data,3)));
                    case 'lin'
                        handles.imHandle = imagesc(prePhasor.data_average);            
                end            
            else
                switch scale
                    case 'log'
                        handle = implay(log10(prePhasor.data));
                    case 'lin'
                        handle = implay(prePhasor.data);
                end
            end
            handle.Visual.ColorMap.MapExpression = 'hot'; 
            handle.Visual.ColorMap.UserRangeMin = 0.1;
            handle.Visual.ColorMap.UserRangeMax = max_val;
            handle.Visual.ColorMap.UserRange = max_val;
        end                       
        
        function handles = show_data_single(prePhasor, scale, index)
            if nargin == 1
                index = 1;
                scale = 'lin';
            elseif nargin == 2
                index = 1;
            end
            handles.figHandle = figure;
            switch scale
                case 'lin'
                    handles.imHandle = imagesc(abs(prePhasor.data(:,:,index)));
                    handles.colorBar = colorbar;
                case 'log'
                    handles.imHandle = imagesc(log10(abs(prePhasor.data(:,:,index))));
                    handles.colorBar = colorbar;
                    ylabel(handles.colorBar,'log');
            end
            axis image;            
            
            colormap jet;
            title([prePhasor.data_meta.sample_name ' | Scan ' num2str(prePhasor.data_meta.scan_number) ' | Frame ' num2str(index)]);
        end                
        
        function handles = show_data_average(prePhasor,scale)
            if ~exist('scale')
                scale = 'lin';
            end
            
            handles.figHandle = figure;            
            
            try 
                switch scale
                    case 'log'
                        handles.imHandle = imagesc(log10(prePhasor.data_average)); 
                        handles.colorBar = colorbar;
                        ylabel(handles.colorBar,'log');
                    case 'lin'
                        handles.imHandle = imagesc(prePhasor.data_average);            
                end
                axis image;            

                colormap jet;
                title(['Average: ' prePhasor.data_meta.sample_name ' | Scan ' num2str(prePhasor.data_meta.scan_number)]);
            catch
                warning('Average data first!');
            end
        end
        
        function handles = show_data_integral(prePhasor)             %should output an appropriate type of a plot
            handles.figHandle = figure;            
            if prePhasor.data_integral
                try
                    handles.imHandle = plot(prePhasor.data_meta.nanomax.gonphi, prePhasor.data_integral,'-o');
                catch
                    handles.imHandle = plot(prePhasor.data_integral,'-o');
                end
                ylabel('Integral intensity');
                xlabel('Scan motor position');
                title([prePhasor.data_meta.sample_name ' | Scan ' num2str(prePhasor.data_meta.scan_number)]);
            elseif ismatrix(prePhasor.data_integral)
                
                % Plotting
                try
                    hVector = (-round(size(prePhasor.data_integral,2)/2):round(size(prePhasor.data_integral,2)/2)-1).*prePhasor.data_meta.nanomax.step_h*1e6;
                    vVector = (round(size(prePhasor.data_integral,1)/2):-1:-(round(size(prePhasor.data_integral,1)/2)-1)).*prePhasor.data_meta.nanomax.step_v*1e6;
                catch
                    hVector = 1:size(prePhasor.data_integral,2);
                    vVector = 1:size(prePhasor.data_integral,1);
                end
                
                try
                    switch scale
                        case 'lin'
                            handles.imHandle = imagesc(hVector,vVector,prePhasor.data_integral);axis image;colormap bone;colorbar;axis xy
                            xlabel('Scan position, [um]');ylabel('Scan position, [um]');
                        case 'log'
                            handles.imHandle = imagesc(hVector,vVector,log10(prePhasor.data_integral));axis image;colormap bone;colorbar;axis xy
                            xlabel('Scan position, [um]');ylabel('Scan position, [um]');
                    end                    
                catch
                    warning('Can not plot an integral map');
                end
                title(['Average: ' prePhasor.data_meta.sample_name ' | Scan ' num2str(prePhasor.data_meta.scan_number)]);
            end
        end   
        
    end
    
    % Save methods
    methods
        function save_gif(prePhasor,user_name)
            disp('Saving GIF animation...');
            file_temp = [prePhasor.data_meta.sample_name,'_',num2str(prePhasor.data_meta.scan_number)];
            mkdir(fullfile(prePhasor.data_meta.save_folder,file_temp));            
            if nargin>1
                gif_name = fullfile(prePhasor.data_meta.save_folder,file_temp,[user_name '.gif']);
            else                
                gif_name = fullfile(prePhasor.data_meta.save_folder,file_temp,[prePhasor.data_meta.sample_name,'_',num2str(prePhasor.data_meta.scan_number),'.gif']);
            end
            f1 = figure;
            for ii = 1:size(prePhasor.data,3)                
                imagesc(log10(squeeze(prePhasor.data(:,:,ii))));
                colormap hot; axis image; %caxis([0.1 0.8]);
                title([prePhasor.data_meta.sample_name ' | Scan ' num2str(prePhasor.data_meta.scan_number) ' | Frame ' num2str(ii)]);           
                GIFanimation(gif_name, f1, 0.1, size(prePhasor.data,3), ii);
            end
            disp('Done!');
        end
        
        function save_mask(prePhasor)
            disp('Saving data to .mat ...');
            file_temp = [prePhasor.data_meta.sample_name,'_',num2str(prePhasor.data_meta.scan_number)];
            mkdir(fullfile(prePhasor.data_meta.save_folder,file_temp));
            
            mask = prePhasor.mask;
            
            if ndims(prePhasor.mask) == 2
                suffix = '2D';
            else
                suffix = '3D';
            end
            
            if nargin>1
                save(fullfile(prePhasor.data_meta.save_folder,file_temp,[user_name,'_mask_',suffix,'.mat']),'mask','-v7.3');
            else                
                save(fullfile(prePhasor.data_meta.save_folder,file_temp,[file_temp,'_mask_',suffix,'.mat']),'mask','-v7.3');            
            end
            disp('Done!');
        end
        
        function save_bin(prePhasor,user_name)           
            disp('Saving data to .bin ...'); 
            file_temp = [prePhasor.data_meta.sample_name,'_',num2str(prePhasor.data_meta.scan_number)];
            mkdir(fullfile(prePhasor.data_meta.save_folder,file_temp));
            data = prePhasor.data;
            if nargin>1
                name = fullfile(prePhasor.data_meta.save_folder,file_temp,[user_name,sprintf('_%d_%d_%d_',size(data,1),size(data,2),size(data,3)),'.bin']);
                fid = fopen(name,'wb');            
            else
                name = fullfile(prePhasor.data_meta.save_folder,file_temp,[file_temp,sprintf('_%d_%d_%d_',size(data,1),size(data,2),size(data,3)),'.bin']);
                fid = fopen(name,'wb');            
            end

            fwrite(fid,data,'double');
            fclose(fid);
            fprintf('Saved: %s \n', name);
        end
        
        function save_data(prePhasor,user_name)
            disp('Saving data to .mat ...');
            file_temp = [prePhasor.data_meta.sample_name,'_',num2str(prePhasor.data_meta.scan_number)];
            mkdir(fullfile(prePhasor.data_meta.save_folder,file_temp));
                        
            if nargin>1
                save(fullfile(prePhasor.data_meta.save_folder,file_temp,[user_name,'.mat']),'scan','-v7.3');
            else
                save(fullfile(prePhasor.data_meta.save_folder,file_temp,[file_temp,'.mat']),'scan','-v7.3');
            end
            
            if prePhasor.crop_flag
                data = prePhasor.data_crop; %#ok<*PROPLC>
            else
                data = prePhasor.data;
            end
            
            if prePhasor.data_meta.save_diff
                if nargin>1
                    save(fullfile(prePhasor.data_meta.save_folder,file_temp,[user_name,'_diff.mat']),'data','-v7.3');
                else                
                    save(fullfile(prePhasor.data_meta.save_folder,file_temp,[file_temp,'_diff.mat']),'data','-v7.3');            
                end
            end
            
            disp('Done!');                        
            
            if strcmp(prePhasor.data_meta.save_formats,'bin')
                      
            end
            
            disp('Done!');
        end
    end        
end

