
% roi_this
% creates ROIs for the current image
% 
% Written by Daniel Buscombe, various times in 2012 and 2013
% while at
% School of Marine Science and Engineering, University of Plymouth, UK
% then
% Grand Canyon Monitoring and Research Center, U.G. Geological Survey, Flagstaff, AZ 
% please contact:
% dbuscombe@usgs.gov
% for lastest code version please visit:
% https://github.com/dbuscombe-usgs
% see also (project blog):
% http://dbuscombe-usgs.github.com/
%====================================
%   This function is part of 'dgs-gui' software
%   This software is in the public domain because it contains materials that originally came 
%   from the United States Geological Survey, an agency of the United States Department of Interior. 
%   For more information, see the official USGS copyright policy at 
%   http://www.usgs.gov/visual-id/credit_usgs.html#copyright
%====================================

if ~isfield(sample(ix),'whole_roi') && sample(ix).whole_roi~=1
    
    % if ~exist('whole_roi','var') || whole_roi
    whole_roi=0;
    sample(ix).num_roi=0;
    for k=1:sample(ix).num_roi
        sample(ix).roi{k}=[];
        sample(ix).roi_x{k}=[];
        sample(ix).roi_y{k}=[];
        sample(ix).roi_line{k}=[];
    end
    chx = get(ax,'Children');
    if length(chx)>=2
        chx(end)=[];
        delete(chx)
    end
    axes(ax)
end


[Nv,Nu,blank] = size(sample(1).data);
set(ax,'xlim',[-2 Nu+2],'ylim',[-2 Nv+2])

sample(ix).whole_roi=0;

sample(ix).num_roi=sample(ix).num_roi+1;

%[blank,rectpos]=crop_image(ax);
rectpos=round(getrect(ax));

% define the points for the line to be drawn
sample(ix).roi_x{sample(ix).num_roi} =round([rectpos(1), rectpos(1)+rectpos(3), rectpos(1)+rectpos(3), ...
    rectpos(1), rectpos(1)]);
sample(ix).roi_y{sample(ix).num_roi} = round([rectpos(2), rectpos(2), rectpos(2)+rectpos(4), ...
    rectpos(2)+rectpos(4), rectpos(2)]);
sample(ix).roi_line{sample(ix).num_roi} = line(sample(ix).roi_x{sample(ix).num_roi},...
    sample(ix).roi_y{sample(ix).num_roi},'color','red','linewidth',2);

sample(ix).roi{sample(ix).num_roi}=sample(ix).data(min(sample(ix).roi_y{sample(ix).num_roi}):...
    max(sample(ix).roi_y{sample(ix).num_roi}),...
    min(sample(ix).roi_x{sample(ix).num_roi}):...
    max(sample(ix).roi_x{sample(ix).num_roi}));

ButtonName = questdlg('Apply this ROI to all images?','ROI', ...
    'Yes','No', 'Yes');

if strcmp(ButtonName,'Yes')
    
    wh = waitbar(0,'Please wait, applying ROIs ...');
    
    for ii=2:length(sample)
        
        % read data in if not already done so
        if isempty(sample(ii).data)
            sample(ii).data=imread([image_path char(image_name(ii))]);
            
            if numel(size(sample(ii).data))==3
                sample(ii).data=double(0.299 * sample(ii).data(:,:,1) + 0.5870 * ...
                    sample(ii).data(:,:,2) + 0.114 * sample(ii).data(:,:,3));
            else
                sample(ii).data=double(sample(ii).data);
            end
            
        end
        im=sample(ii).data;
        [n,m,p] = size(im);
        
        
        try
            [n,m,p] = size(im);
            
            v = ver;
            if any(strcmp('Statistics Toolbox', {v.Name}))
                % cosine taper
                w = .25;
                window = repmat(tukeywin(n,w),1,m).*rot90(repmat(tukeywin(m,w),1,n));
                
                for i = 1:p
                    im(:,:,i) = im(:,:,i).*window;
                end
            end
        catch
            continue
        end
        
        sample(ii).data=im;
        
        
        % first remove previous rois
        if sample(ii).num_roi>0 % && sample(ii).whole_roi~=1
            for k=1:sample(ii).num_roi
                sample(ii).roi{k}=[];
                sample(ii).roi_x{k}=[];
                sample(ii).roi_y{k}=[];
                sample(ii).roi_line{k}=[];
            end
        end
        
        try
            
            sample(ii).num_roi=1;
            
            % define the points for the line to be drawn
            sample(ii).roi_x{sample(ii).num_roi} =round([rectpos(1), rectpos(1)+rectpos(3), rectpos(1)+rectpos(3), ...
                rectpos(1), rectpos(1)]);
            sample(ii).roi_y{sample(ii).num_roi} = round([rectpos(2), rectpos(2), rectpos(2)+rectpos(4), ...
                rectpos(2)+rectpos(4), rectpos(2)]);
            sample(ii).roi_line{sample(ii).num_roi} = line(sample(ii).roi_x{sample(ii).num_roi},...
                sample(ii).roi_y{sample(ii).num_roi},'color','red','linewidth',2);
            
            sample(ii).roi{sample(ii).num_roi}=sample(ii).data(min(sample(ii).roi_y{sample(ii).num_roi}):...
                max(sample(ii).roi_y{sample(ii).num_roi}),...
                min(sample(ii).roi_x{sample(ii).num_roi}):...
                max(sample(ii).roi_x{sample(ii).num_roi}));
        catch
            disp(['Image ',num2str(ii),' is a different size, ROI not applied'])
            for k=1:sample(ii).num_roi
                sample(ii).roi{k}=[];
                sample(ii).roi_x{k}=[];
                sample(ii).roi_y{k}=[];
                sample(ii).roi_line{k}=[];
                sample(ii).num_roi=0;
            end
        end
        waitbar(ii/length(sample),wh)
        
    end
    close(wh)
    
    
end

% for ii=1:length(sample)
%     sample(ii).roi(cellfun(@isempty,sample(ii).roi))=[];
%     sample(ii).roi_x(cellfun(@isempty,sample(ii).roi_x))=[];
%     sample(ii).roi_y(cellfun(@isempty,sample(ii).roi_y))=[];
%     sample(ii).roi_line(cellfun(@isempty,sample(ii).roi_line))=[];    
% end


set(findobj('tag','current_image'),'userdata',sample);

if ~isempty(sample(ix).dist)
    uiwait(msgbox('... remember to calculate size distribution again!','New ROIs defined ...','modal'));
end

clear rectpos
