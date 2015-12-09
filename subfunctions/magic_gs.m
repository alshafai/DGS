
function [s1, s2, s3]=magic_gs(auto)

R = 0.5;

% this next section fits an ellipse to the contour R in auto
k=size(auto,1);
ind=[-(k/2)+1:k/2]; % centre portion
[xind,yind]=find(auto>R); % find data greater than R

try
P=[ind(xind)',ind(yind)']; % To reduce the computation time, work with the boundary points only
catch
f=find(yind>max(ind)); yind(f)=[]; xind(f)=[];
P=[ind(xind)',ind(yind)']; % To reduce the computation time, work with the boundary points only
end

try
K = convhulln(P);
K = unique(K(:));
PK = P(K,:)';

[d N] = size(PK);
Q = zeros(d+1,N);
Q(1:d,:) = PK(1:d,1:N);
Q(d+1,:) = ones(1,N);

% initializations
count = 1;
err = 1;
u = (1/N) * ones(N,1);          % 1st iteration
tolerance=.01;

% use the Khachiyan Algorithm (see Ref 3)
while err > tolerance,
    X = Q * diag(u) * Q';       % X = \sum_i ( u_i * q_i * q_i')  is a (d+1)x(d+1) matrix
    M = diag(Q' * inv(X) * Q);  % M the diagonal vector of an NxN matrix
    [maximum j] = max(M);
    step_size = (maximum - d -1)/((d+1)*(maximum-1));
    new_u = (1 - step_size)*u ;
    new_u(j) = new_u(j) + step_size;
    count = count + 1;
    err = norm(new_u - u);
    u = new_u;
end

% (x-c)' * A * (x-c) = 1
% It computes a dxd matrix 'A' and a d dimensional vector 'c' as the center
% of the ellipse. 
U = diag(u);

% the A matrix for the ellipse
A = (1/d) * inv(PK * U * PK' - (PK * u)*(PK*u)' );
% matrix contains all the information regarding the shape of the ellipsoid

% center of the ellipse 
c = PK * u;

[U Q V] = svd(A);

r1 = 1/sqrt(Q(1,1));
r2 = 1/sqrt(Q(2,2));
v = [r1 r2 c(1) c(2) V(1,1)]';

% get ellipse points
N = 100; % number of points on ellipse
dx = 2*pi/N; % spacing
theta = v(5); %orinetation
Rot = [ [ cos(theta) sin(theta)]', [-sin(theta) cos(theta)]']; % rotation matrix
Xe=zeros(1,N); Ye=zeros(1,N); % pre-allocate
   
for i = 1:N
     ang = i*dx;
     x = v(1)*cos(ang);
     y = v(2)*sin(ang);
     d1 = Rot*[x y]';
     Xe(i) = d1(1) + v(3);
     Ye(i) = d1(2) + v(4);
end   

% mean euclidean distance
h=sum(imag(sqrt((0-Xe.^2)+(0-Ye.^2))))/length(Xe);

s1=(2*pi)/(v(1)^-1); 
s2=(2*pi)/(v(2)^-1);
s3=(2*pi)/(h^-1);

catch
    a = rad_av(auto);
    s1 = (2*pi)*find(a<R, 1,'first');
    s2 = s1; s3 = s1;
end
