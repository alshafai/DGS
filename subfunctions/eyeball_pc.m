
% eyeball_pc
function C=eyeball_pc(im,n,split)
close all

warning off all
% function to carry out 'point counts' of grain diameters on an image of sediment
%
% INPUTS
% im = digital image or image filename
% n= number of grains you wish to count
% split - either 1, 4 or 16. The number of equal-sized windows the image will
% be split into in order to zoom into the grains
% choose 1 if the grains are large compared to image area, 16 if the grains
% are very small, and 4 if somewhere in the middle
%
% for the counts, the choice of protocol is left open:
% to randomise objectively, two options are given:
% a grid is drawn
% green dots, the coordinates from a random number generator, are plotted
%
% in order to identify grains, the figure window is automatically maximised
%
% grains are identified by the user using a cursor. red lines are drawn across the grains 
%
% OUTPUT
% C= vector of counted grain sizes, in pixels
%
%  Author:  Daniel Buscombe
%           School of Marine Science & Engineering
%           University of Plymouth, Drake Circus, Plymouth, Devon PL48AA, UK
%           daniel.buscombe@plymouth.ac.uk
%  Version: Beta        Revision: 16 April, 2010

if ~isnumeric(im)
im=double(rgb2gray(imread(deblank(im))));
end

if nargin==1 % i.e. if just the image is given as input
    n=100; % deafaults to 100 grains
    split=4; % 4 sections
end


if split==1
    
        C=[];
    
        imagesc(im), axis image
        colormap gray
        grid on
        set(gca,'GridLineStyle','-')
        hold on, 
        maximise
        
        %pnt=rand(n,2);
        %plot(pnt(:,1).*size(p,2),pnt(:,2).*size(p,1),'g.','MarkerSize',16)
        for i=1:n
            a=ginput(2); % two points
            C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
            line(a(:,1),a(:,2),'color','r','LineWidth',2)
        end

        close

elseif split ==4
    
    C=[];
    e=round(n/4);
    ex=n-(e*4);
    
    [m,n,o]=size(im);
    if o>1
        im=col2gray(im);
    end
    
    if m~=n
    s=min(m,n);
    if rem(s,2)==1
        s=s-1; % if odd number
    end
    im=im(1:s,1:s);   % if image is not square, resize to 
    % [1:smallest even dimension]
    end
    
    p=size(im,1)/2; % split image into 4 quarters 
    q=size(im,1)/2;
    [ m, n ] = size(im);
    if rem(m, p) || rem(n, q)
      error(sprintf([ 'Can''t create %d-by-%d submatrices from' ...
                      ' a %d-by-%d matrix.' ], p, q, m, n));
    end

    y1 = reshape(im, [ p m/p q n/q ]);
    y1 = permute(y1, [ 1 3 2 4 ]); 
    l1=size(y1);
    
    for ind=1:2
        for jnd=1:2
        
        p=y1(:,:,ind,jnd);
 
        imagesc(p), axis image
        colormap gray
        grid on
        set(gca,'GridLineStyle','-')
        hold on, 
        maximise
        
        if ind==2 && jnd==2
        pnt=rand(e+ex,2);
        plot(pnt(:,1).*size(p,2),pnt(:,2).*size(p,1),'g.','MarkerSize',16)
        for i=1:e+ex
            a=ginput(2);
            C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
            line(a(:,1),a(:,2),'color','r','LineWidth',2)
        end
        else
        pnt=rand(e,2);
        plot(pnt(:,1).*size(p,2),pnt(:,2).*size(p,1),'g.','MarkerSize',16)
        for i=1:e
            a=ginput(2);
            C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
            line(a(:,1),a(:,2),'color','r','LineWidth',2)
        end
        end 

        close
        
        end
    end


    elseif split ==16

e=round(n/16);
ex=n-(e*16);
disp(['warning - you will do ',num2str(ex),' extra points in the first window'])

[m,n]=size(im);

sec=im(1:round(m)/4,1:round(n)/4);
C=[];
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 1 out of 16')
maximise

p=rand(e+ex,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e+ex
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 2
sec=im(round(m)/4:round(m)/2,1:round(m)/4);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 2 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 3
sec=im(round(m)/2:[round(m)/4]*3,1:round(m)/4);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 3 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 4
sec=im([round(m)/4]*3:m,1:round(m)/4);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 4 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close


%---------------- 5
sec=im(1:round(m)/4,round(n)/4:round(n)/2);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 5 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 6
sec=im(round(m)/4:round(m)/2,round(m)/4:round(n)/2);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 6 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 7
sec=im(round(m)/2:[round(m)/4]*3,round(m)/4:round(n)/2);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 7 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 8
sec=im([round(m)/4]*3:m,round(m)/4:round(n)/2);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 8 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close


%---------------- 9
sec=im(1:round(m)/4,round(n)/2:(round(n)/4)*3);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 9 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 10
sec=im(round(m)/4:round(m)/2,round(n)/2:(round(n)/4)*3);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 10 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 11
sec=im(round(m)/2:[round(m)/4]*3,round(n)/2:(round(n)/4)*3);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 11 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 12
sec=im([round(m)/4]*3:m,round(n)/2:(round(n)/4)*3);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 12 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close



%---------------- 13
sec=im(1:round(m)/4,3*(round(n)/4):n);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 13 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 14
sec=im(round(m)/4:round(m)/2,3*(round(n)/4):n);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 14 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 15
sec=im(round(m)/2:[round(m)/4]*3,3*(round(n)/4):n);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 15 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

%---------------- 16
sec=im([round(m)/4]*3:m,3*(round(n)/4):n);
imagesc(sec), axis image
colormap gray
grid on
set(gca,'GridLineStyle','-')
hold on, title('Section 16 out of 16')
maximise

p=rand(e,2);
plot(p(:,1).*size(sec,2),p(:,2).*size(sec,1),'y.','MarkerSize',16)

for i=1:e
    a=ginput(2);
    C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
    line(a(:,1),a(:,2),'color','r','LineWidth',2)
end
close

else

            C=[];
    
        imagesc(im), axis image
        colormap gray
        grid on
        set(gca,'GridLineStyle','-')
        hold on, 
        maximise
        
        pnt=rand(n,2);
        plot(pnt(:,1).*size(p,2),pnt(:,2).*size(p,1),'g.','MarkerSize',16)
        for i=1:n
            a=ginput(2);
            C=[C;ceil(sqrt(diff(a(:,1))^2+diff(a(:,2))^2))];
            line(a(:,1),a(:,2),'color','r','LineWidth',2)
        end

        close
    
end

C=ceil(C(:))+(10.*(ceil(C(:))./100)); % add 10% to grain sizes 
%(natural bias is to count egde to edge, rather than pore to pore)

warning on all


function maximise

units=get(gcf,'units');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gcf,'units',units)



