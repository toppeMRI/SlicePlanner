function pfile2hdf(pfile, echo, readoutfile, hdffile)
% function pfile2hdf(pfile, [echo=1, readoutfile='readout.mod', hdffile='Localizer.h5])
%
% Recon a 3D spin-warp (GRE) Pfile acquired with TOPPE, 
% and save (magnitude) image stack in HDF5 format that can be loaded into the GUI in ../GUI/.

if ~exist('hdffile', 'var')
	hdffile = 'Localizer.h5';
end

if ~exist('readoutfile', 'var')
	readoutfile = 'readout.mod';
end

if ~exist('echo', 'var')
	echo = 1;
end

if strcmp(pfile, 'test')
	sub_test();
	return
end

% zero-fill factor
zf = 2;                  

% voxel dimensions
seq.fov = 24;  % cm (isotropic)
seq.n = 120;   % matrix size (isotropic)
dx = seq.fov/seq.n/zf;      % cm
dy = dx;
dz = dx;

% get Pfile header
[~, rdb_hdr] = toppe.utils.loadpfile(pfile);

% reconstruct (magnitude) images
[~, imsos] = toppe.utils.recon3dft(pfile, 'echo', echo, 'readoutFile', readoutfile, ...
    'alignWithUCS', true);

% make matrix large to increase image on screen in Java GUI
dat = fftshift(fftn(fftshift(imsos)));
[nx,ny,nz] = size(dat);
dzf = zeros(zf*nx, zf*ny, zf*nz);
dzf((end/2-nx/2+1):(end/2+nx/2),(end/2-ny/2+1):(end/2+ny/2),(end/2-nz/2+1):(end/2+nz/2)) = dat;
imstmp = abs(fftshift(ifftn(fftshift(dzf))));

% pad image with zeros to make matrix size isotropic (makes it easier to program the GUI)
[nx,ny,nz] = size(imstmp);
nmax = max([ny ny nz]);
imsos = zeros(nmax, nmax, nmax);
imsos((end/2-nx/2+1):(end/2+nx/2),(end/2-ny/2+1):(end/2+ny/2),(end/2-nz/2+1):(end/2+nz/2)) = imstmp;

% get new size
[nx,ny,nz] = size(imsos);

% magnitude image
imsos = abs(imsos);

% make brightness across image more uniform
%imsos = imsos/max(imsos(:));
%imsos = imsos.*(1-0.9*imsos);

% scale to [0,255] Integer range (will be converted to RGB in Java GUI)
imsos = round(imsos/max(imsos(:))*255);

% write to HDF5 file
% Images are written as 2D slices, in axial/sagittal/coronal orientations
% HDF5 likes to transpose matrices, so transpose each 2D image before writing to file.
system(sprintf('rm -f %s', hdffile));

hdf5write(hdffile, '/Dims/nx', int16(nx));
hdf5write(hdffile, '/Dims/ny', int16(ny), 'writemode', 'append');
hdf5write(hdffile, '/Dims/nz', int16(nz), 'writemode', 'append');

hdf5write(hdffile, '/Voxel/dx', dx, 'writemode', 'append');
hdf5write(hdffile, '/Voxel/dy', dy, 'writemode', 'append');
hdf5write(hdffile, '/Voxel/dz', dz, 'writemode', 'append');

imsos = uint8(round(imsos));

imsosAxi = imsos;  
for iz = 1:nz
	hdf5write(hdffile, sprintf('/Ax/slice%d', iz), imsos(:,:,iz), 'writemode', 'append');
end

imsosSag = permute(imsos, [2 3 1]);
for ix = 1:nx
	hdf5write(hdffile, sprintf('/Sag/slice%d', ix), imsosSag(:,:,ix), 'writemode', 'append');
end

imsosCor = permute(imsos, [1 3 2]);
for iy = 1:ny
	hdf5write(hdffile, sprintf('/Cor/slice%d', iy), imsosCor(:,:,iy), 'writemode', 'append');
end

% write each FID (view) separately to avoid confusion!
if 0
hdffile = 'LocScanFIDs.h5';
hdf5write(hdffile, '/Dims/nx', int16(nx));
hdf5write(hdffile, '/Dims/ny', int16(ny), 'writemode', 'append');
hdf5write(hdffile, '/Dims/nz', int16(nz), 'writemode', 'append');
for iy = 1:ny
	for iz = 1:nz
		hdf5write(hdffile, sprintf('/FID/y%d/z%d', iy, iz), imsos(:,iy,iz), 'writemode', 'append');
	end
end
end

% write 3D volume to file
%hdf5write('LocScan3D.h5', '/Dims', int16([nx ny nz]));
%hdf5write('LocScan3D.h5', '/Image', uint8(imsos), 'writemode', 'append');

return


function sub_test()

datdir = '/export/data/jfnielse/stack-of-spirals-presto-bold-fmri/fbirn,14Aug2019/';

datdir = '/export/data/jfnielse/stack-of-spirals-presto-bold-fmri/invivo,fieldmap,30Aug2019/';
pfile = [datdir 'P,fmri,b0,30Aug2019.7'];

datdir = '/export/data/jfnielse/stack-of-spirals-presto-bold-fmri/gephantom,30Aug2019/';
pfile = [datdir 'P,gephantom,30Aug2019,b0.7'];

readoutfile = [datdir 'readout_withheader.mod'];

hdffile = 'Localizer.h5';
echo = 1;
pfile2hdf(pfile, echo, readoutfile, hdffile);

data = h5read(hdffile, '/Sag/slice36');
im(data);

return
