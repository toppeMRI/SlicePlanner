# 3D SPGR localizer (scout) sequence for SlicePlanner


## Requirements

The TOPPE Matlab toolbox

```
$ git clone git@github.com:toppeMRI/toppe.git
```

In Matlab:
```
>> addpath('./toppe/');
```

## Create the TOPPE sequence files

In Matlab:
```
>> makelocscan;       % creates localizer.tar
```

## Execute the TOPPE sequence on a GE scanner

1. Untar the file 'localizer.tar' in /usr/g/research/pulseq/loc/ on the scanner host.
2. Prescribe the toppev4 interpreter sequence. Set the number of slices to 120 or greater. Freq. Dir: R/L.
3. Start scan.

## Convert P-file to 'Localizer.h5'

Reconstruct a 3D image volume and write it to an hdf5 file:

```
>> pfile = 'P12345.7';
>> readoutFile = 'readout.mod';
>> pfile2hdf(pfile, 1, readoutFile);  
```



