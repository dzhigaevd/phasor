classdef Phasor < handle
    %   PHASOR Summary of this class goes here
    %   3D version
    % h - horizontal axis of the diffraction pattern
    % v - vertical axis of the diffraction pattern
    % r - rocking curve direction
    
    properties
        data;
        reciprocal_coordinates;
        real_coordinates;
        data_meta;
        experiment;
        data_binned; 
        data_gap;
        support;
        object;
        object_fft_mod        
        probe;        
        metric;
        snr;
           
        % Constants. They are needed for correct labeling of axes
        h                       = 4.1357e-15;                                  % Plank's constant
        c                       = 2.99792458e8;                                % Speed of light in vacuum
    end
    
    
    % - Constructor of the object -
    
    methods
        function obj = Phasor(input_param)                
            obj.data_meta = input_param;
            obj.metric.val.sharpness = [];
            obj.metric.val.reciprocal = [];
            
            obj.experiment.wavelength = obj.h*obj.c/obj.data_meta.energy;
                        
            disp('PHASOR object created succesfully');
        end
    end  
    
    methods
        % Load methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function load_mat(obj)
            % Dynamic determination of the dataset name
            in = load(obj.data_meta.data_path); 
            names = fieldnames(in);
            temp = in.(names{1});
            clear in;           
            
            temp(isnan(temp)) = 0;
            if temp<0
                obj.data_gap = (temp==-1);
                temp(obj.data_gap) = 0;
                disp('Gaps are recorded!');
            else
                obj.data_gap = (temp==-1);
                disp('No gaps');
            end            
            
            obj.data = single(sqrt(double(temp)));
            obj.data_meta.data_size = size(obj.data);
            
            fprintf('Data loaded: %s\n',obj.data_meta.data_path);
            fprintf('Data size: [%d %d %d]\n',obj.data_meta.data_size);
        end                              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % General functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function ram2gpu(obj)
            gpuDevice(1);
            obj.data   = gpuArray(obj.data);
            obj.support = gpuArray(obj.support);
            obj.object = gpuArray(obj.object);   
            obj.metric.val.reciprocal = gpuArray(obj.metric.val.reciprocal);
            obj.metric.val.sharpness = gpuArray(obj.metric.val.sharpness);
            disp('Calculating with GPU')
        end
        
        function gpu2ram(obj)
            obj.object  = gather(obj.object);
            obj.support = gather(obj.support);
            obj.data    = gather(obj.data);
            obj.metric.val.reciprocal = gather(obj.metric.val.reciprocal);
            obj.metric.val.sharpness = gather(obj.metric.val.sharpness);
            disp('Arrays gathered from GPU')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        % Class methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function create_support(obj,path)
            % Creat initial support            
            if nargin == 1
                % Creation of the support first
                autoCorr = abs(fftNc(obj.data.^2));   
                obj.support = (autoCorr./max(autoCorr(:))) > 0.1;
                disp('Auto-correlation support created!')
            elseif nargin == 2  
                try
                    obj.support = load(path);
                    disp('External support loaded!\n');
                catch
                    error('Cannot load the support! Check the path!');
                end
            end
        end
        
        function create_object(obj,type)
            if isempty(obj.support)
                error('No support defined! Create support first!')
            else
                if nargin == 1
                    type = 'flat';
                end

                % Initialize the starting guess
                % options: 
                % amp  - random amplitude, constant phase
                % ph   - random phase, constant amplitude
                % amph - random phase, random amplitude
                % flat - constant amplitude
                if strcmp(type,'amp')
                    obj.object = rand(obj.data_meta.data_size).*obj.support;
                elseif strcmp(type,'ph')
                    obj.object = ones(obj.data_meta.data_size).*exp(1j*(-pi+2*pi.*rand(obj.data_meta.data_size))).*obj.support;
                elseif strcmp(type,'amph')
                    obj.object = rand(obj.data_meta.data_size).*exp(1j*(-pi+2*pi.*rand(obj.data_meta.data_size))).*obj.support;
                elseif strcmp(type,'flat')
                    obj.object = ones(obj.data_meta.data_size).*obj.support;
                end

                fprintf('Object of type %s created!\n',type);
            end
        end
        
        function treshold_data(obj,threshold)            
            obj.data(obj.data<threshold) = 0;
        end
        
        function out = PS(obj)
            % Support projection operator
            out = obj.object.*obj.support;% g(k+1) = g'(k) for support region  
        end

        function out = PM(obj)
            % Modulus projection operator
%             if obj.data_meta.rac % apply Regularized Amplitude Constraint
%                 inFT = fftNc(param.O);
%                 Bufferdata = inFT(param.Bmstop==1);
%                 switch param.RACtype
%                     case 'UNI'
%                         E_l = 1-((1./param.data).^2)/24;                
%                     case 'GAUSS'
%                         E_l = 1-((1./param.data).^2)/8;
%                     case 'POISS'
%                         E_l = 1-((1./param.data).^1.199)*0.236;
%                 end
%                 E_l(isinf(E_l)) = 1;
%                 inFT = E_l.*param.data.*param.dataMask.*inFT./(abs(inFT));
%                 inFT(param.Bmstop==1) = Bufferdata;
%                 out = ifftNc(inFT); % g'(k)                
%             else
                inFT = fftNc(obj.object);
                Bufferdata = inFT(obj.data_gap==1);
                inFT = obj.data.*inFT./(abs(inFT));
                inFT(obj.data_gap==1) = Bufferdata;
                out = ifftNc(inFT); % g'(k)        
%             end                    
        end
                
        function calculate_metric(obj)
            if isempty(obj.metric.val)
                val_pos = 1;
            else
                switch obj.data_meta.metric.type
                    case 'reciprocal' 
                        val_pos = numel(obj.metric.val.reciprocal)+1;
                    case 'sharpness'
                        val_pos = numel(obj.metric.val.sharpness)+1;
                    case 'both'
                        val_pos = numel(obj.metric.val.sharpness)+1;                        
                end
            end
            
            switch obj.data_meta.metric.type
                case 'reciprocal'                                        
                    obj.metric.val.reciprocal(1,val_pos)  = sum(sum(sum((abs(fftNc(obj.object.*obj.support)).*abs(obj.data_gap-1)-obj.data).^2)))/...
                        sum(sum(sum((obj.data).^2)));             
                case 'real'                    
                    obj.metric.val(1,val_pos) = sum(sum(sum(abs(obj.object.*abs(obj.support-1)).^2)))/...
                        sum(sum(sum(abs(obj.object.*obj.support).^2))); 
                    % Real space error for HIO
                case 'sharpness'                    
                    d = (abs(obj.object)).^4;
                    obj.metric.val.sharpness(1,val_pos) = sum(d(:));   
                case 'both'
                    obj.metric.val.reciprocal(1,val_pos)  = sum(sum(sum((obj.object_fft_mod.*abs(obj.data_gap-1)-obj.data).^2)))/...
                        sum(sum(sum((obj.data.^2))));
                    d = (abs(obj.object)).^4;
                    obj.metric.val.sharpness(1,val_pos) = sum(d(:)); 
                    
            end            
        end
        
        function calculate_snr(obj)
            try
                % from "Invariant error metrics for image reconstruction" J. R. Fienup
                obj.snr = (1-mean(obj.metric.val.reciprocal))/mean(obj.metric.val.reciprocal);
                fprintf('Approxiamte SNR: %.3f\n',obj.snr);
            catch
                error('Chi squared has to be calculated!')
            end
        end
        
        function calculate_object_fft_mod(obj)
            % Should we compare the same volume inside the support
                % only, since everything outside is not an object?
            obj.object_fft_mod = abs(fftNc(obj.object.*obj.support));
        end
             
        % Phase retrieval methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function SW(obj, sigma, treshold)
            if isempty(obj.object)
                error('No object defined! Create object first!')
            else
                if numel(obj.data_meta.data_size) == 3
                    temp = imgaussfilt3(abs(obj.object),[sigma sigma sigma]);
                elseif numel(obj.data_meta.data_size) == 2
                    temp = imgaussfilt(abs(obj.object),[sigma sigma]);
                else
                    warning('Wrong dimensionality of the input data');
                end
                                
                obj.support = temp>treshold*max(temp(:)); 
                fprintf('Support is shrinked by sigma %.3f at treshold %.3f\n',sigma, treshold);
            end
        end
        
        function ER(obj,n,show_progress)
            for ii = 1:n                                                                
                obj.object         = PM(obj); % Magnitude projection operator
                obj.object         = obj.object.*obj.support; % Support projection operator   
                
                obj.calculate_object_fft_mod;
                obj.calculate_metric;
                
                if show_progress     
                    fprintf('ER iteration: %d\n',ii);
                end
            end             
        end
        
        function HIO(obj,n,beta,show_progress)
            if nargin == 2
                beta = 0.9;
                show_progress = 0;
            elseif nargin == 3
                show_progress = 0;
            end            
            
            for ii = 1:n                                             
                buffer = PM(obj); % Modulus projection operator                                
                obj.object = buffer.*obj.support+(1-obj.support).*(obj.object-beta*buffer);                                 
                
                obj.calculate_object_fft_mod;
                obj.calculate_metric;
                
                if show_progress 
                    fprintf('HIO iteration: %d\n',ii);
                end
            end             
        end
        
        function SF(obj,n,show_progress)            
            for ii = 1:n   
                
                obj.object = (2*obj.support-1).*PM(obj); % Magnitude projection operator  
                
                obj.calculate_object_fft_mod;
                obj.calculate_metric;
                
                if show_progress
                    fprintf('SF iteration: %d\n',ii); 
                end
            end             
        end
        
        function RAAR(obj,n,beta,show_progress)            
            for ii = 1:n                           
                buffer = PM(obj);
                obj.object = 0.5*beta*( (2*obj.support-1).*(2*buffer-obj.object)+obj.object)+(1-beta)*buffer; % Magnitude projection operator   
                
                obj.calculate_object_fft_mod;
                obj.calculate_metric;
                
                if show_progress   
                    fprintf('RAAR iteration: %d\n',ii); 
                end
            end             
        end
        
        % Post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function remove_phase_ramp(obj)
            c0 = size(obj.object)/2;                                                   
            if numel(size(obj.object)) == 3 
                for ii = 1:3
                    dummy = fftNc(obj.object);
                    c1 = ndimCOM(abs(dummy).^4,'auto');                
                    d =((-1)^ii)*(c1-c0);                                       
                    fprintf('Data center: [%.2f, %.2f, %.2f]\nShift by [%.2f, %.2f, %.2f]\n', c1(1), c1(2), c1(3),d(1),d(2),d(3));                    
                    dummy = imtranslate(dummy,d);
                    dummy = ifftNc(dummy);
                    % Average phase value
                    obj.object = dummy.*exp(-1j*mean(angle(dummy(:))));
                end
            elseif numel(obj.data_meta.data_size) == 2
                c1 = ndimCOM(obj.data,'auto');
                fprintf('Data center: [%.2f, %.2f]', c1(1), c1(2));
            else
                warning('Data dimensionality is wrong!');
            end                    
            disp('Phase ramp removed!');
        end
        
        function center_data(obj)
            obj.data = center_array_com(obj.data);
            disp('Data centered!');
        end
        
        function transform2lab(obj,beamline)
            switch beamline
                case '34DC'
                    
                case 'P10' 
            end
        end
        
        function add_pi(obj)
            obj.object = obj.object.*exp(1j*pi);
            disp('+pi value is added to the phase of the object')
        end
        
        % Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_metric(obj)                                                
            if strcmp(obj.data_meta.metric.type,'reciprocal')
                figure;
                ax = semilogy(obj.metric.val.reciprocal,'-k')   ;         
                title(sprintf('Minimum value: %.3f',min(obj.metric.val.reciprocal(:))));
                ylabel('Reciprocal difference'); ax.YGrid = 'on';
            elseif strcmp(obj.data_meta.metric.type,'sharpness')
                figure;
                ax = semilogy(obj.metric.val.sharpness,'-k')  ;          
                title(sprintf('Maximum value: %.3f',max(obj.metric.val.sharpness(:))));             
                ylabel('Sharpness'); ax.YGrid = 'on';
            elseif strcmp(obj.data_meta.metric.type,'both')
                figure;
                subplot(2,1,1);
                ax1 = semilogy(obj.metric.val.reciprocal,'-k');           
                title(sprintf('Minimum value: %.5f',min(obj.metric.val.reciprocal(:))));
                ylabel('Reciprocal difference'); grid on;
                xlabel('Iteration number');
                
                subplot(2,1,2);
                ax2 = semilogy(obj.metric.val.sharpness,'-k');           
                title(sprintf('Maximum value: %.3e',max(obj.metric.val.sharpness(:))));             
                ylabel('Sharpness'); grid on;                
            end
            xlabel('Iteration number');
        end
        
        function slice3d_object(obj)
%             vis3d(double(angle(obj.object).*obj.support, obj.support),'jet','Phase inside support');
            vis3d(double(angle(obj.object).*obj.support), obj.support);
        end
                
        function iso3d_object(obj,input,cmap)
            % Use an input parameter to show other complex valued matrix
            if nargin == 1
                cmap = 'jet';
                input = abs(obj.object)./max(max(max(abs(obj.object))));       
%                 input = input.*.*obj.support;
            elseif nargin == 2
                cmap = 'jet';                
            end            
            handle = figure;            
            panel = uipanel('Parent',handle,'Title','Amp_Phase','FontSize',...
            12,'Units','Normalized','Position',[.1 0.1 .77 .85]);         
            ax = axes('Parent',panel);
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.5,'Units','Normalized',...
            'Position', [0.1 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReal); 
            isoVal = 0.5;
            drawIsosurface(input,isoVal,cmap);
            
            function drawIsosurface(input,isoVal,cmap)
                cla(ax);
                axes(ax);
                isosurface(input,isoVal,angle(obj.object));
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
                rotate3d on;
                grid on;
                axis tight;
                axis equal; 
                axis vis3d;
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];           
                colormap(cmap);
%                 col = colorbar('Parent',handle);
%                 col.Label.String = 'Phase, [rad]';
            end
            
            function slideIsosurfaceReal(hObj,callbackdata)
                isoVal = get(hObj,'Value');                 
                drawIsosurface(input,isoVal,cmap);
            end  
        end
        
        function iso3d_data(obj,input_flag,cmap)
            % Use an input parameter to show fft mod of the object
            if nargin == 1
                cmap = 'jet';
                input = log10(obj.data)./max(log10(obj.data(:)));
                title_val = 'Data';
            elseif nargin == 2
                cmap = 'jet';
                input = log10(obj.object_fft_mod)./max(log10(obj.object_fft_mod(:)));
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
                axis vis3d;
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];            
                colormap(cmap);                 
            end
            
            function slideIsosurfaceReal(hObj,callbackdata)
                isoVal = get(hObj,'Value');                 
                drawIsosurface(input,isoVal);
            end  
        end
        % Save functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function save(obj,marker)
            mkdir('reconstructions');
            save(sprintf('reconstructions/reconstruction_%s.mat',marker),'obj');
        end
    end
end
