

function a = rad_av(auto)
% returns a 1D radial average of 2D autocorrelogram
 
[m,n] = size(auto);
[X Y] = meshgrid(-m/2:m/2-1, -m/2:m/2-1);% Make Cartesian grid
[theta rho] = cart2pol(X, Y); % Convert to polar coordinate axes
rho = round(rho);
i = cell(floor(m/2) + 1, 1);
for r = 0:floor(m/2)
    i{r + 1} = find(rho == r);
end
% radial average of autocorrelation
a = zeros(1, floor(m/2)+1);
warning off
for r = 0:floor(m/2)
    a(1, r + 1) = mean(auto( i{r+1} ) );
end
a(isnan(a)) = []; % remove trailing NaN
warning on