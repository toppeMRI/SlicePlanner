# 3D SPGR 'localizer' sequence for SlicePlanner


## Requirements

The TOPPE Matlab toolbox

```
$ git clone git@github.com:toppeMRI/toppe.git
```

In Matlab:
```
>> addpath('./toppe/');
```

## Create 3D spin-warp sequence in TOPPE format

In Matlab:
```
>> main;       % creates toppev3,localizer.tgz
```

## Execute the sequence on a GE scanner

1. Untar the file 'toppev3,localizer.tgz' in /usr/g/bin/ on the scanner host.
2. Prescribe the toppev3 interpreter sequence. Set the number of slices to 120 or greater. 
3. Start scan.

## Convert P-file to 'Localizer.h5'

```
>> pfile = 'P12345.7';
>> readoutFile = 'readout.mod';
>> pfile2hdf(pfile, 1, readoutFile);  
```



