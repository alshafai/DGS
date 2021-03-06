
% dgs_gui_swopsimages
% callback for main program, swops images
% 
% Written by Daniel Buscombe, various times in 2012 - 2014
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

% get current image
ix=get(findobj('tag','PickImage'),'value');
% and current image data
sample=get(findobj('tag','current_image'),'userdata');

% % get next image
if isempty(sample(ix).data)
    sample(ix).data=imread([image_path char(image_name(ix))]);
    
    if numel(size(sample(ix).data))==3
        sample(ix).data=double(0.299 * sample(ix).data(:,:,1) + 0.5870 * ...
            sample(ix).data(:,:,2) + 0.114 * sample(ix).data(:,:,3));
    else
        sample(ix).data=double(sample(ix).data);
    end
    
    im=sample(ix).data;
    
    try
        [n,m,p] = size(im);
        
%         v = ver;
%         if any(strcmp('Statistics Toolbox', {v.Name}))
%             % cosine taper
%             w = .25;
%             window = repmat(tukeywin(n,w),1,m).*rot90(repmat(tukeywin(m,w),1,n));
%             
%             for i = 1:p
%                 im(:,:,i) = im(:,:,i).*window;
%             end
%         end
    catch
        continue
    end
    
    %     [n,m,p] = size(im);
    %     % cosine taper
    %     w = .25; % width of cosine in percent of width of X
    %     window = repmat(tukeywin(n,w),1,m).*rot90(repmat(tukeywin(m,w),1,n));
    %
    %     for i = 1:p
    %         im(:,:,i) = im(:,:,i).*window;
    %     end
    sample(ix).data=im;
    
end

set(findobj('tag','current_image'),'cdata',sample(ix).data);

% set resolution bar to whatever the resolution actually is
set(findobj('tag','res'),'String',num2str(sample(ix).resolution));

% remove plotted roi lines from previous image
chx = get(ax,'Children');
if length(chx)>=2
    chx(end)=[];
    delete(chx)
end

h=findobj('tag','current_image');
set(h,'cdata',sample(ix).data); % make first image appear

[Nv,Nu,blank] = size(sample(ix).data);
set(h,'xdata',1:Nu); % scales and labels
set(h,'ydata',1:Nv);

set(ax,'ylim',[1,size(sample(ix).data,1)])
set(ax,'xlim',[1,size(sample(ix).data,2)])

% if navigating back, draw roi lines back on
if sample(ix).num_roi>0
    for k=1:sample(ix).num_roi
        sample(ix).roi_line{k} = line(sample(ix).roi_x{k},sample(ix).roi_y{k},'color','red','linewidth',2);
    end
end

% first set axes ticks to be increments of 500
set(ax,'ytick',linspace(1,size(sample(ix).data,1),2))
set(ax,'xtick',linspace(1,size(sample(ix).data,2),2))
% scale current x and y labels
set(ax,'xticklabels',num2str(get(ax,'xtick')'.*sample(ix).resolution))
set(ax,'yticklabels',num2str(get(ax,'ytick')'.*sample(ix).resolution))
% axis tight

if isfield(sample(ix),'roi_line')
    if iscell(sample(ix).roi_line)
        if ~isempty(sample(ii).roi_line)
            sample(ix).roi_line{1} = line(sample(ix).roi_x{1},...
                sample(ix).roi_y{1},'color','red','linewidth',5);
        end
    end
end

if isempty(sample(ix).auto)
    
    chx = get(ax3,'Children');
    if length(chx)>=2
        chx(end)=[];
        delete(chx)
    end
    axes(ax3)
    title('')
    
    [Nv,Nu,blank] = size(sample(ix).data);
    
    % calculate 2D autocorrel
    im=sample(ix).data(1:min(Nu,Nv),1:min(Nu,Nv));
    % 2D-FFT transform on de-meaned image
    % power spectrum
    mag=abs(fft2(fftshift(im-mean(im(:))))).^2;
    %Shift zero-frequency component to centre of spectrum
    auto=fftshift(real(ifft2(mag)));
    auto = auto./max(auto(:));
    
    [centx,centy] = find(auto==1);
    % spectify number of lags to compute
    l = length(auto);
    nlags=round(l/8);
    % centre 2d autocorrelogram
    auto = auto(centx-nlags:centx+nlags,centy-nlags:centy+nlags);
    
    sample(ix).auto = auto;
    [Nv,blank,blank] = size(sample(ix).auto);
    
    h=findobj('tag','auto_image');
    
    set(h,'userdata',sample);
    set(h,'cdata',sample(ix).auto); % make fi
    
    set(findobj('tag','auto_axes'),'xlim',[-2 2+Nv],...
        'ylim',[-2 2+Nv])
    grid off
    title('2D autocorrelation')
    axes(ax)
    
end
%         cla(ax3)
if ~isempty(sample(ix).dist)
    
    %     cla(ax3)
    chx = get(ax3,'Children');
    if length(chx)>=2
        chx(end)=[];
        delete(chx)
    end
    
    h=findobj('tag','auto_image');
    
    tmpimage=sample(ix).data; %roi{1};
    [Nv,Nu,blank] = size(tmpimage);
    tmpimage=tmpimage(round((Nv/2)-sample(ix).percentiles(8)*1/sample(ix).resolution):...
        round((Nv/2)+sample(ix).percentiles(8)*1/sample(ix).resolution),...
        round((Nu/2)-sample(ix).percentiles(8)*1/sample(ix).resolution):...
        round((Nu/2)+sample(ix).percentiles(8)*1/sample(ix).resolution));
    [Nv,Nu,blank] = size(tmpimage);
    set(h,'cdata',tmpimage); % make fi
    axes(ax3)
    set(findobj('tag','auto_axes'),'xlim',[-2 2+Nv],...
        'ylim',[-2 2+Nv])
    set(ax3,'xticklabels',num2str(get(ax3,'xtick')'.*sample(ix).resolution))
    set(ax3,'yticklabels',num2str(get(ax3,'ytick')'.*sample(ix).resolution))
    
    grid off
    title('Sample Of Image')
    hold on
    plot([Nv/2 Nv/2],...
        [Nu/2 (Nu/2)+sample(ix).percentiles(5)*1/sample(ix).resolution],'r-','linewidth',2)
    text(Nv/2,Nu/2,'D_{50}','color','g','fontsize',12)
    
else
    chx = get(ax3,'Children');
    if length(chx)>=2
        chx(end)=[];
        delete(chx)
    end
    [Nv,blank,blank] = size(sample(ix).auto);
    
    h=findobj('tag','auto_image');
    
    set(h,'userdata',sample);
    set(h,'cdata',sample(ix).auto); % make fi
    
    set(findobj('tag','auto_axes'),'xlim',[-2 2+Nv],...
        'ylim',[-2 2+Nv])
    grid off
    axes(ax3)
    title('2D autocorrelation')
end
axes(ax)


if ~isempty(sample(ix).dist)
    
    cla(ax2)
    
    h=findobj('Tag','plot_axes');
    axes(h)
    bar(sample(ix).dist(:,1),sample(ix).dist(:,2))
    if sample(ix).resolution==1
        xlabel('Size (Pixels)')
    else
        xlabel('Size (mm)')
    end
    ylabel('Density')
    axis tight
    set(gca,'ydir','normal')
    text(.7,.92,['Mean = ',num2str(sample(ix).arith_moments(1),3)],'units','normalized','fontsize',7)
    text(.7,.85,['Sorting = ',num2str(sample(ix).arith_moments(2),3)],'units','normalized','fontsize',7)
    text(.7,.78,['Skewness = ',num2str(sample(ix).arith_moments(3),3)],'units','normalized','fontsize',7)
    text(.7,.70,['D_{10} = ',num2str(sample(ix).percentiles(2),3)],'units','normalized','fontsize',7)
    text(.7,.62,['D_{50} = ',num2str(sample(ix).percentiles(5),3)],'units','normalized','fontsize',7)
    text(.7,.54,['D_{90} = ',num2str(sample(ix).percentiles(8),3)],'units','normalized','fontsize',7)
    grid off
    axes(ax2)
    title('Size Distribution')
    axes(ax)
    
else
    cla(ax2)
    title('')
end


set(findobj('tag','current_image'),'userdata',sample);

clear chx k n Nu Nv mag im auto nlags l centx centy tmpimage h


% update title
% a=get(findobj('tag','im_axes1'),'title');
% set(get(findobj('tag','im_axes1'),'title'),'string',char(sample(ix).name));


