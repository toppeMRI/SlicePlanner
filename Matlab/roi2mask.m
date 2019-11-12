function mask = roi2mask(roi,seq)
%
% Input: 
%  roi    struct    Created with this GUI. Load with toppe.getroi()
%                   Distance (spatial) unit is pixels, for now. 
%  seq    struct    Sequence parameters. Get with ../getParams()
%
% Output:
%  mask    [nx ny nz]   Binary mask corresponding to input ROI

[nx,ny,nz] = deal(seq.n, seq.n, seq.n);

mask = zeros(nx,ny,nz);

zf = seq.zeroFillFactor;

% width, height, thickness
w = round(roi.w/zf);
h = round(roi.h/zf);
t = round(roi.t/zf);

% center location (offset from iso-center)
x = roi.x/zf;    % readout (1st Matlab matrix dimension)
y = roi.y/zf;    % y phase-encode direction
z = roi.z/zf;    % z phase-encode ("slice") direction

% rotation matrix
rotmat = roi.rotmat;

% loop over pixels inside ROI
for ii = -w/2:w/2
	for jj = -h/2:h/2
		for kk = -t/2:t/2
			% p = vector representing location of one pixel
			p = roi.rotmat*[ii jj kk]';        % rotate 
			p = p + [x y z]';                  % translate 
			p = round(p + [nx/2 ny/2 nz/2]');  % convert to matrix index

			% update mask. Include neghboring pixels to avoid a 'perforated' mask
			%mask(p(2),p(1),p(3)) = 1;   % NB! Note x/y transpose, since GUI +x = 2nd matrix dimension in Matlab
			for ni = -1:1
				for nj = -1:1
					for nk = -1:1
						ptmp = p + [ni nj nk]';

						% don't go beyond the 3D matrix
						ptmp = max(ptmp,1);
						ptmp = min(ptmp,seq.n);
						mask(ptmp(1),ptmp(2),ptmp(3)) = 1;
					end
				end
			end
		end
	end
end

return

