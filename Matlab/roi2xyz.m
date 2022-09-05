function [X, Y, Z] = roi2xyz(roi, FOV, showROI)
% function [X, Y, Z] = roi2xyz(roi, [FOV, showROI])
%
% Get spatial locations X, Y, Z (in cm) for grid locations
% inside a rectangular ROI (that is possibly off-center and rotated).
%
% Inputs:
%   roi   struct containing rectangular ROI parameters. See toppe.getroi()
%   FOV   [1 3]   field of view (cm). Temporary fix since roi is currently in pixel units
%
% Output:
%   X/Y/Z   [nVoxInt 1]   Voxels locations (cm) in the interior of the ROI.

if nargin < 2
    FOV = [24 24 20];  % cm
end
if nargin < 3
    showROI = true;
end

N = [roi.wmax roi.tmax roi.hmax];  % matrix size

% Distance between center of edge points
% (max voxel distance from center of FOV is FOVc/2)
FOVc = FOV - FOV./N;

% At the moment, roi units are in pixels (TODO: make roi units = cm).
% So here's a fix for now.
roi.w = roi.w * FOV(1)/N(1);  % width (cm) (x)
roi.h = roi.h * FOV(2)/N(2);  % height  (z)
roi.t = roi.t * FOV(3)/N(3);  % thickness  (y)
roi.x = roi.x * FOV(1)/N(1);  % center (x coordinate)
roi.y = roi.y * FOV(2)/N(2);  % center
roi.z = roi.z * FOV(3)/N(3);  % center

[X, Y, Z] = ndgrid(linspace(-1,1,N(1))*FOVc(1)/2, ...
                   linspace(-1,1,N(2))*FOVc(2)/2, ...
                   linspace(-1,1,N(3))*FOVc(3)/2);

% get voxel locations inside rectangle (before rotating and translating)
mask = ones(N);
mask(X > roi.w/2 | X < -roi.w/2) = 0;
mask(Y > roi.t/2 | Y < -roi.t/2) = 0;
mask(Z > roi.h/2 | Z < -roi.h/2) = 0;
X = X(logical(mask));  % column vector
Y = Y(logical(mask));
Z = Z(logical(mask));

% rotate
LOCS = roi.rotmat * [X Y Z]';  % roi.rotmat = 3x3 rotation matrix
X = LOCS(1,:);
Y = LOCS(2,:);
Z = LOCS(3,:);

% translate
X = X + roi.x;
Y = Y + roi.y;
Z = Z + roi.z;

% remove locations outside FOV
inds = X > FOVc(1)/2 | X < -FOVc(1)/2 ...
     | Y > FOVc(2)/2 | Y < -FOVc(2)/2 ...
     | Z > FOVc(3)/2 | Z < -FOVc(3)/2;
X(inds) = [];
Y(inds) = [];
Z(inds) = [];

% display
if showROI
    % convert voxel locations (in cm) to matrix indeces
    Xinds = N(1)/2 +  round(X/FOV(1)*N(1));     
    Yinds = N(2)/2 +  round(Y/FOV(2)*N(2));    
    Zinds = N(3)/2 +  round(Z/FOV(3)*N(3));    

    % create mask and display
    mask = zeros(N);
    for ii = 1:length(Xinds)
        mask(Xinds(ii), Yinds(ii), Zinds(ii)) = 1;
    end

    im(mask);   % from MIRT
    title('NB! Display may contain zeros (holes) due to rounding');
end



