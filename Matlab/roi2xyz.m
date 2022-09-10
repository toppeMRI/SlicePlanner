function [X, Y, Z] = roi2xyz(roi, N, FOV, showROI)
% function [X, Y, Z] = roi2xyz(roi, N, FOV, [showROI=true])
%
% Get spatial locations X, Y, Z (in cm) for grid locations
% inside a rectangular ROI (that is possibly off-center and rotated).
% Grid locations are defined by N and FOV.
%
% Inputs:
%   roi     struct containing rectangular ROI parameters. See toppe.getroi()
%   N       [1 3]   matrix size
%   FOV     [1 3]   field of view (cm)
%   shoROI  true/false   
%
% Output:
%   X/Y/Z   [nVoxInt 1]   Voxels locations (cm) in the interior of the ROI.

% internal calculations are in mm
FOV = FOV*10;   % mm

if nargin < 4
    showROI = true;
end

% Distance between center of edge points
% (max voxel distance from center of FOV is FOVc/2)
FOVc = FOV - FOV./N;

[X, Y, Z] = ndgrid(linspace(1,-1,N(1))*FOVc(1)/2, ...
                   linspace(1,-1,N(2))*FOVc(2)/2, ...
                   linspace(1,-1,N(3))*FOVc(3)/2);

% get voxel locations inside rectangle (before rotating and translating)
mask = ones(N);
mask(X > roi.w/2 | X < -roi.w/2) = 0;
mask(Y > roi.h/2 | Y < -roi.h/2) = 0;
mask(Z > roi.t/2 | Z < -roi.t/2) = 0;
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
% this is approximate since some grid points may be missed after rotating and translating
if showROI
    % convert voxel locations (in mm) to matrix indeces
    % Note negative sign due to universal coordinate system 
    Xinds = N(1)/2 +  round(-X/FOV(1)*N(1));     
    Yinds = N(2)/2 +  round(-Y/FOV(2)*N(2));    
    Zinds = N(3)/2 +  round(-Z/FOV(3)*N(3));    

    % create mask and display
    mask = zeros(N);
    for ii = 1:length(Xinds)
        mask(Xinds(ii), Yinds(ii), Zinds(ii)) = 1;
    end

    im(mask);   % from MIRT
    title('NB! Display may contain zeros (holes) due to rounding');
end

% convert to cm
X = X(:)/10;
Y = Y(:)/10;
Z = Z(:)/10;

