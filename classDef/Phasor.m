classdef Phasor < handle
    %   PHASOR Summary of this class goes here
    %   3D version
    % h - horizontal axis of the diffraction pattern
    % v - vertical axis of the diffraction pattern
    % r - rocking curve direction
    
    properties
        data;
        dataTime;
        reciprocal_coordinates;
        real_coordinates;
        data_meta;
        experiment;
        data_binned; 
        data_gap;
        support;
        object;
        mcf; % mutual coherence function
        object_fft_mod;
        objectFT_buffer;
        reconstruction;
        probe;        
        metric;
        snr;
        current_iterate;
        
        % Constants. They are needed for correct labeling of axes
        h                       = 4.1357e-15;                                  % Plank's constant
        c                       = 2.99792458e8;                                % Speed of light in vacuum
    end    
    
    % - Constructor of the object -    
    methods
        function phasor = Phasor(input_param,data)       
            phasor.data_meta = input_param;
            phasor.metric.val.sharpness = [];
            phasor.metric.val.reciprocal = [];
            phasor.data_meta.algo_list = '';
            phasor.experiment.wavelength = phasor.h*phasor.c/phasor.data_meta.energy;
            phasor.current_iterate = 1;
            disp('PHASOR object created succesfully!');
            if nargin==2                                
                phasor.data = single(sqrt(double(data)));
                phasor.data_meta.data_size = size(phasor.data);
            
                fprintf('Data loaded: %s\n',phasor.data_meta.data_path);
                fprintf('Data size: [%d %d %d]\n',phasor.data_meta.data_size);      
                
                phasor.data(isnan(phasor.data)) = 0;
                if phasor.data<0
                    phasor.data_gap = (phasor.data==-1);
                    phasor.data(phasor.data_gap) = 0;
                    disp('Gaps are recorded!');
                else
                    phasor.data_gap = (phasor.data==-1);
                    disp('No gaps');
                end  
            end
            
        end
    end  
    
    methods
        % Load methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function load_mat(phasor)
            % Dynamic determination of the dataset name
            in = load(phasor.data_meta.data_path); 
            names = fieldnames(in);
            temp = in.(names{1});
            clear in;           
            
            temp(isnan(temp)) = 0;
            if temp<0
                phasor.data_gap = (temp==-1);
                temp(phasor.data_gap) = 0;
                disp('Gaps are recorded!');
            else
                phasor.data_gap = (temp==-1);
                disp('No gaps');
            end            
            
            phasor.data = single(sqrt(double(temp)));
            phasor.data_meta.data_size = size(phasor.data);
            
            fprintf('Data loaded: %s\n',phasor.data_meta.data_path);
            fprintf('Data size: [%d %d %d]\n',phasor.data_meta.data_size);
        end                              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % General functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        function ram2gpu(phasor)
            gpuDevice(1);
            phasor.data   = gpuArray(phasor.data);
            phasor.support = gpuArray(phasor.support);
            phasor.object = gpuArray(phasor.object);   
            phasor.object_fft_mod = gpuArray(phasor.object_fft_mod);
            phasor.metric.val.reciprocal = gpuArray(phasor.metric.val.reciprocal);
            phasor.metric.val.sharpness = gpuArray(phasor.metric.val.sharpness);
            disp('Calculating with GPU')
        end
        
        function gpu2ram(phasor)
            phasor.object  = gather(phasor.object);
            phasor.object_fft_mod = gather(phasor.object_fft_mod);
            phasor.support = gather(phasor.support);
            phasor.data    = gather(phasor.data);
            phasor.metric.val.reciprocal = gather(phasor.metric.val.reciprocal);
            phasor.metric.val.sharpness = gather(phasor.metric.val.sharpness);
            disp('Arrays gathered from GPU')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        % Class methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function create_support(phasor,path)
            % Creat initial support            
            if nargin == 1
                % Creation of the support first
                autoCorr = abs(fftNc(phasor.data.^2));   
                phasor.support = (autoCorr./max(autoCorr(:))) > 0.1;
                disp('Auto-correlation support created!')
            elseif nargin == 2  
                try
                    phasor.support = load(path);
                    disp('External support loaded!\n');
                catch
                    error('Cannot load the support! Check the path!');
                end
            end
        end
        
        function create_object(phasor,type)
            if isempty(phasor.support)
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
                    phasor.object = rand(phasor.data_meta.data_size).*phasor.support;
                elseif strcmp(type,'ph')
                    phasor.object = ones(phasor.data_meta.data_size).*exp(1j*(-pi+2*pi.*rand(phasor.data_meta.data_size))).*phasor.support;
                elseif strcmp(type,'amph')
                    phasor.object = rand(phasor.data_meta.data_size).*exp(1j*(-pi+2*pi.*rand(phasor.data_meta.data_size))).*phasor.support;
                elseif strcmp(type,'flat')
                    phasor.object = ones(phasor.data_meta.data_size).*phasor.support;
                end
                
                if phasor.data_meta.partial_coherence
%                    phasor.mcf = (11,11,11,0.5,0.5,0.5,1);
%                    phasor.mcf = arr
                end
                fprintf('Object of type %s created!\n',type);
            end
        end
        
        function threshold_data(phasor,threshold)            
            phasor.data(phasor.data<threshold) = 0;
        end
        
        function out = PS(phasor)
            % Support projection operator
            out = phasor.object.*phasor.support;% g(k+1) = g'(k) for support region  
        end

        function out = PM(phasor)
            % Modulus projection operator
%             if phasor.data_meta.rac % apply Regularized Amplitude Constraint
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
            if ~phasor.data_meta.partial_coherence                 
                objectFT = fftNc(phasor.object);                
                Bufferdata = objectFT(phasor.data_gap==1);
                objectFT = phasor.data.*objectFT./(abs(objectFT));
                objectFT(phasor.data_gap==1) = Bufferdata;
                out = ifftNc(objectFT); % g'(k)
            elseif phasor.data_meta.partial_coherence...
                    && (phasor.current_iterate <= phasor.data_meta.partial_coherence_start...
                    || phasor.current_iterate >= phasor.data_meta.partial_coherence_end)
                objectFT = fftNc(phasor.object);
                if phasor.current_iterate == phasor.data_meta.partial_coherence_start-1
                    phasor.objectFT_buffer = objectFT;
                end                                
                Bufferdata = objectFT(phasor.data_gap==1);
                objectFT = phasor.data.*objectFT./(abs(objectFT));
                objectFT(phasor.data_gap==1) = Bufferdata;
                out = ifftNc(objectFT); % g'(k)                
            elseif phasor.data_meta.partial_coherence...
                    && phasor.current_iterate >= phasor.data_meta.partial_coherence_start...
                    && phasor.current_iterate <= phasor.data_meta.partial_coherence_end
                
                objectFT = fftNc(phasor.object);
                % Update the mcf                                
                IdeltaK = 2*abs(objectFT).^2-abs(phasor.objectFT_buffer).^2;
                IdeltaK(IdeltaK <0) =0;
                
                for ii = 1:phasor.data_meta.pc_iterations
                    fprintf('Update of coherence function %d\n',ii);
                    mcf = phasor.mcf;
                    phasor.mcf = mcf*convn(flip(flip(flip(IdeltaK,1),2),3),phasor.data.^2./(convn(IdeltaK,mcf)));
                end
                
                phasor.objectFT_buffer = objectFT;
                % Create the current estimate of  of the Ipc
                Ipc = convn(abs(objectFT).^2,phasor.mcf);
                objectFT = phasor.data.*objectFT./sqrt(Ipc);
                out = ifftNc(objectFT); % g'(k)
            end                    
        end
                
        function calculate_metric(phasor)
            if isempty(phasor.metric.val)
                val_pos = 1;
            else
                switch phasor.data_meta.metric.type
                    case 'reciprocal' 
                        val_pos = numel(phasor.metric.val.reciprocal)+1;
                    case 'sharpness'
                        val_pos = numel(phasor.metric.val.sharpness)+1;
                    case 'both'
                        val_pos = numel(phasor.metric.val.sharpness)+1;                        
                end
            end
            
            switch phasor.data_meta.metric.type
                case 'reciprocal'                                        
                    phasor.metric.val.reciprocal(1,val_pos)  = sum(sum(sum((abs(fftNc(phasor.object.*phasor.support)).*abs(phasor.data_gap-1)-phasor.data).^2)))/...
                        sum(sum(sum((phasor.data).^2)));             
                case 'real'                    
                    phasor.metric.val(1,val_pos) = sum(sum(sum(abs(phasor.object.*abs(phasor.support-1)).^2)))/...
                        sum(sum(sum(abs(phasor.object.*phasor.support).^2))); 
                    % Real space error for HIO
                case 'sharpness'                    
                    d = (abs(phasor.object)).^4;
                    phasor.metric.val.sharpness(1,val_pos) = sum(d(:));   
                case 'both'
                    phasor.metric.val.reciprocal(1,val_pos)  = sum(sum(sum((phasor.object_fft_mod.*abs(phasor.data_gap-1)-phasor.data).^2)))/...
                        sum(sum(sum((phasor.data.^2))));
                    d = (abs(phasor.object)).^4;
                    phasor.metric.val.sharpness(1,val_pos) = sum(d(:)); 
                    
            end            
        end
        
        function calculate_snr(phasor)
            try
                % from "Invariant error metrics for image reconstruction" J. R. Fienup
                phasor.snr = (1-mean(phasor.metric.val.reciprocal))/mean(phasor.metric.val.reciprocal);
                fprintf('Approxiamte SNR: %.3f\n',phasor.snr);
            catch
                error('Chi squared has to be calculated!')
            end
        end
        
        function calculate_object_fft_mod(phasor,input)
            % Should we compare the same volume inside the support
                % only, since everything outside is not an object?
            if nargin == 1
                phasor.object_fft_mod = abs(fftNc(phasor.object.*phasor.support));
            else
                phasor.object_fft_mod = abs(fftNc(input.*phasor.support));
            end
        end

        function truncate_amplitude(phasor, threshold)
            phasor.support = abs(phasor.object)>threshold*(max(abs(phasor.object(:))));
        end

        % Phase retrieval methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function SW(phasor, sigma, threshold)
            if isempty(phasor.object)
                error('No object defined! Create object first!')
            else
                if numel(phasor.data_meta.data_size) == 3
                    temp = imgaussfilt3(abs(phasor.object),[sigma sigma sigma]);
                elseif numel(phasor.data_meta.data_size) == 2
                    temp = imgaussfilt(abs(phasor.object),[sigma sigma]);
                else
                    warning('Wrong dimensionality of the input data');
                end
                                
                phasor.support = temp>threshold*max(temp(:)); 
                
                phasor.calculate_object_fft_mod;
                phasor.calculate_metric;                
                fprintf('Support is shrinked by sigma %.3f at threshold %.3f\n',sigma, threshold);
            end
        end
        
        function final_SW(phasor, sigma, threshold_range, steps)
            if isempty(phasor.object)
                error('No object defined! Create object first!')
            else
                if numel(phasor.data_meta.data_size) == 3
                    temp = imgaussfilt3(abs(phasor.object),[sigma sigma sigma]);
                elseif numel(phasor.data_meta.data_size) == 2
                    temp = imgaussfilt(abs(phasor.object),[sigma sigma]);
                else
                    warning('Wrong dimensionality of the input data');
                end
                
                threshold = linspace(threshold_range(1),threshold_range(2),steps);
                
                for ii = 1:numel(threshold)
                    phasor.support = temp>threshold(ii)*max(temp(:));
                    for kk = 1:5
                        buffer         = PM(phasor); % Magnitude projection operator
                        buffer         = buffer.*phasor.support; % Support projection operator 
                    end
                    
                    phasor.calculate_object_fft_mod(buffer);
                    phasor.calculate_metric;
                    err_val(ii) = phasor.metric.val.reciprocal(end);
                end                
                [val, pos] = find(err_val == min(err_val));
                figure; plot(err_val);
                phasor.support = temp>threshold(pos(1))*max(temp(:));
                phasor.calculate_object_fft_mod;
                phasor.calculate_metric;
                
                fprintf('Support is shrinked by sigma %.3f at threshold %.3f\n',sigma, threshold(pos(1)));
            end
        end
        

        
        function ER(phasor,n,show_progress)
            for ii = 1:n                                                                
                phasor.object         = PM(phasor); % Magnitude projection operator
                phasor.object         = phasor.object.*phasor.support; % Support projection operator   
                
                phasor.calculate_object_fft_mod;
                phasor.calculate_metric;
                
                phasor.print_iteration;
                phasor.current_iterate = phasor.current_iterate+1;
            end
            
            if show_progress     
                fprintf('+ %d ER iterations\n',ii);
            end
                
            if ~ismember('ER',phasor.data_meta.algo_list)
                phasor.data_meta.algo_list = [phasor.data_meta.algo_list,'ER_'];
            end
        end
        
        function HIO(phasor,n,HIO_beta,show_progress)
            if nargin == 2
                HIO_beta = 0.9;
                show_progress = 0;
            elseif nargin == 3
                show_progress = 0;
            end            
            
            for ii = 1:n          
                
                buffer = PM(phasor); % Modulus projection operator                                
                phasor.object = buffer.*phasor.support+(1-phasor.support).*(phasor.object-HIO_beta*buffer);                                 
                
                phasor.calculate_object_fft_mod;
                phasor.calculate_metric;
                
                phasor.print_iteration;
                phasor.current_iterate = phasor.current_iterate+1;
            end 
            
            if show_progress     
                fprintf('+ %d HIO iterations\n',ii);
            end
            
            if ~ismember('HIO',phasor.data_meta.algo_list)
                phasor.data_meta.algo_list = [phasor.data_meta.algo_list,sprintf('HIO_',n)];
            end
        end
        
        function SF(phasor,n,show_progress)            
            for ii = 1:n   
                
                phasor.object = (2*phasor.support-1).*PM(phasor); % Magnitude projection operator  
                
                phasor.calculate_object_fft_mod;
                phasor.calculate_metric;  
                
                phasor.print_iteration;
                phasor.current_iterate = phasor.current_iterate+1;
            end   
            
            if show_progress     
                fprintf('+ %d SF iterations\n',ii);
            end
            
            if ~ismember('SF',phasor.data_meta.algo_list)
                phasor.data_meta.algo_list = [phasor.data_meta.algo_list,'SF'];
            end
        end
        
        function RAAR(phasor,n,RAAR_beta,show_progress)            
            for ii = 1:n                           
                buffer = PM(phasor);
                phasor.object = 0.5*RAAR_beta*( (2*phasor.support-1).*(2*buffer-phasor.object)+phasor.object)+(1-RAAR_beta)*buffer; % Magnitude projection operator   
                
                phasor.calculate_object_fft_mod;
                phasor.calculate_metric;
                
                phasor.current_iterate = phasor.current_iterate+1;
            end   
            
            if show_progress     
                fprintf('+ %d RAAR iterations\n',ii);
            end
            
            if ~ismember('RAAR',phasor.data_meta.algo_list)
                phasor.data_meta.algo_list = [phasor.data_meta.algo_list,'RAAR_'];
            end
        end
        
        % Post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function remove_phase_ramp(phasor)
            c0 = size(phasor.object)/2;                                                   
            if numel(size(phasor.object)) == 3 
                for ii = 1:3
                    dummy = fftNc(phasor.object);
                    c1 = ndimCOM(abs(dummy).^4,'auto');                
                    d =((-1)^ii)*(c1-c0);                                       
                    fprintf('Data center: [%.2f, %.2f, %.2f]\nShift by [%.2f, %.2f, %.2f]\n', c1(1), c1(2), c1(3),d(1),d(2),d(3));                    
                    dummy = imtranslate(dummy,d);
                    dummy = ifftNc(dummy);
                    % Average phase value
                    phasor.object = dummy.*exp(-1j*mean(angle(dummy(:))));
                end
            elseif numel(phasor.data_meta.data_size) == 2
                c1 = ndimCOM(phasor.data,'auto');
                fprintf('Data center: [%.2f, %.2f]', c1(1), c1(2));
            else
                warning('Data dimensionality is wrong!');
            end                    
            disp('Phase ramp removed!');
        end
        
        function center_data(phasor)
            [phasor.data,center] = center_array_com(phasor.data);
            fprintf('Data centered at [%.2f %.2f %.2f]\n',(center));
        end
        
        function center_object(phasor)
            [phasor.object, com_val] = center_array_com(phasor.object);
            phasor.support = center_array_com(phasor.support, com_val);
            disp('Object centered!');
        end
        
        function transform2lab(phasor)
            switch phasor.data_meta.beamline
                case '34idc'
                    DCS_to_SS;
                case 'P10' 
            end
        end
        
        function add_pi(phasor)
            phasor.object = phasor.object.*exp(1j*pi);
            disp('+pi value is added to the phase of the object')
        end
        
        % Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function print_iteration(phasor)
            clc;
            fprintf('### Current iteration %d ###\n', phasor.current_iterate);
        end
        
        function plot_metric(phasor)                                                
            if strcmp(phasor.data_meta.metric.type,'reciprocal')
                figure;
                ax = semilogy(phasor.metric.val.reciprocal,'-k')   ;         
                title(sprintf('Minimum value: %.3f',min(phasor.metric.val.reciprocal(:))));
                ylabel('Reciprocal difference'); ax.YGrid = 'on';
            elseif strcmp(phasor.data_meta.metric.type,'sharpness')
                figure;
                ax = semilogy(phasor.metric.val.sharpness,'-k')  ;          
                title(sprintf('Maximum value: %.3f',max(phasor.metric.val.sharpness(:))));             
                ylabel('Sharpness'); ax.YGrid = 'on';
            elseif strcmp(phasor.data_meta.metric.type,'both')
                figure;
                subplot(2,1,1);
                ax1 = semilogy(phasor.metric.val.reciprocal,'-k');           
                title(sprintf('Minimum value: %.5f',min(phasor.metric.val.reciprocal(:))));
                ylabel('Reciprocal difference'); grid on;
                xlabel('Iteration number');
                
                subplot(2,1,2);
                ax2 = semilogy(phasor.metric.val.sharpness,'-k');           
                title(sprintf('Maximum value: %.3e',max(phasor.metric.val.sharpness(:))));             
                ylabel('Sharpness'); grid on;                
            end
            xlabel('Iteration number');
        end
        
        function slice3d_object(phasor)
            vis3d(double(angle(phasor.object).*phasor.support), phasor.support);
        end
        
        function iso3d(phasor)
            % Use an input parameter to show other complex valued matrix            
            cmap = 'jet';
            realObject = phasor.object./max(max(max(abs(phasor.object))));       
            realObject = realObject.*phasor.support;
            
            handle = figure;            
            sH1 = subplot(1,2,1);                             
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.5,'Units','Normalized',...
            'Position', [0.1 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReci); 
            isoVal = 0.5;
            drawIsosurfaceReci(phasor.data./max(phasor.data(:)),isoVal,cmap,sH1);
            
            sH2 = subplot(1,2,2);                          
            uicontrol('Parent',handle,'Style',...
            'slider','Min',0,'Max',1,...
            'Value',0.5,'Units','Normalized',...
            'Position', [0.6 0.05 0.3 0.03],...
            'Callback', @slideIsosurfaceReal); 
            isoVal = 0.5;
            drawIsosurfaceReal(realObject,isoVal,cmap,sH2);

            linkprop([sH1, sH2],{ 'CameraPosition'});
            
            function drawIsosurfaceReal(input,isoVal,cmap,ax)
                cla(ax);
                axes(ax);
                isosurface(abs(input),isoVal,angle(input));
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
                rotate3d on;
                grid on;
                axis tight;
                axis equal; 
                axis vis3d;
                h3 = light; h3.Position = [-1 -1 -1];  
                h4 = light; h4.Position= [1 1 1];           
                colormap(cmap);
            end
            
            function drawIsosurfaceReci(input,isoVal,cmap,ax)
                cla(ax);
                axes(ax);
                isosurface(input,isoVal*max(input(:)));
                xlabel('q_x'); ylabel('q_y'); zlabel('q_z'); 
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
                drawIsosurfaceReal(realObject,isoVal,cmap,sH2);
            end  
            
            function slideIsosurfaceReci (hObj,callbackdata)
                isoVal = get(hObj,'Value');                 
                drawIsosurfaceReci(phasor.data,isoVal,cmap,sH1);
            end  
        end
        
        function iso3d_object(phasor,input,cmap)
            % Use an input parameter to show other complex valued matrix
            if nargin == 1
                cmap = 'jet';
                input = phasor.object./max(max(max(abs(phasor.object))));       
                input = input.*phasor.support;
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
                isosurface(abs(input),isoVal,angle(input));
                xlabel('x, [nm]'); ylabel('y, [nm]'); zlabel('z, [nm]'); 
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
                drawIsosurface(input,isoVal,cmap);
            end  
        end
        
        function iso3d_data(phasor,input_flag,cmap)
            % Use an input parameter to show fft mod of the object
            if nargin == 1
                cmap = 'jet';
                input = (phasor.data)./max((phasor.data(:)));
                title_val = 'Data';
            elseif nargin == 2
                cmap = 'jet';
                input = (phasor.object_fft_mod)./max((phasor.object_fft_mod(:)));
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
        
        % Save functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save the current instance of the object
        function save_object(phasor,marker)            
            rec_path = [phasor.data_meta.save_path, sprintf('/Scan_%s_%s',phasor.data_meta.algo_list,phasor.data_meta.dataTime)];
            mkdir(rec_path);
            save([rec_path, sprintf('/phasor_%s.mat', marker)],'phasor');
            fprintf('Saved the object to: %s',[rec_path, sprintf('/phasor_%s.mat',marker)]);
        end
        
        % save the reconstruction in the current state
        function save_reconstruction(phasor,marker)            
            rec_path = [phasor.data_meta.save_path, sprintf('/Scan_%s_%s',phasor.data_meta.algo_list,phasor.data_meta.dataTime)];
            mkdir(rec_path);
            obj = phasor.object;            
            save([rec_path, sprintf('/reconstruction_%s.mat', marker)],'obj');
            fprintf('Saved the object to: %s',[rec_path, sprintf('/reconstruction_%s.mat\n',marker)]);
        end               
        
        % save vtk file for advanced visualisation in ParaView
        function save2vtk(phasor,marker)
            if ~isempty(phasor.reconstruction)
                rec_path = [phasor.data_meta.save_path, sprintf('/Scan_%s_%s',phasor.data_meta.algo_list,phasor.data_meta.dataTime)];
                mkdir(rec_path);
                savemat2vtk([rec_path, sprintf('/reconstruction_%s.vtk',marker)],phasor.reconstruction.object,phasor.reconstruction.object_pitch);
            else
                warning('Reconstruction is not transformed into sample coordinate system!\n')
                rec_path = [phasor.data_meta.save_path, sprintf('/Scan_%s_%s',phasor.data_meta.algo_list,phasor.data_meta.dataTime)];
                mkdir(rec_path);
                savemat2vtk([rec_path, sprintf('/reconstruction_%s.vtk',marker)],phasor.object);
            end
        end
    end
end

