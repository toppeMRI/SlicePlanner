# SliceRx

## Installing

See [Install.md](Install.md)


## Using

1. Create localizer image volume (or request a sample image from Jon)
	1. Create TOPPE localizer scan files
		```
		>> cd ./experiment/Localizer
		>> main;
		```
	1. Scan :)
	1. Convert the Pfile to an HDF5 file (named 'Localizer.h5') containing the 3D image volume:
		```
		>> cd ../LocalizerScan
		>> pfile2hdf('P12345.7');
		```
1. Start GUI
	```
	$ cd ../GUI
	$ ln -s ../LocalizerScan/Localizer.h5 .
	$ ./start
	```
1. Adjust ROIs, then press 'Export ROIs' which will create the file 'ROI.h5'
1. To load an existing ROI.h5 file, press 'Load ROIs'
1. Load an ROI into Matlab and use it in your pulse design code:
	```
	>> roiId = 1;                                                % which of the ROIs to load
	>> roi = toppe.getroi('ROI.h5', roiId);
	>> rotmat = roi.rotmat;                                      % 3x3 rotation matrix
	>> resLoc = 0.1;                                             % Voxel size of localizer volume (assumed to be isotropic) (cm).            
	>> sliceOffset = resLoc*roi.scanPlaneToIsocenterDistance;    % cm
	```

