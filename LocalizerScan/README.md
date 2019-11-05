# 3D SPGR 'localizer' sequence for toppev3

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
>> main;
```

## Execute the sequence on a GE scanner

Prescribe the toppev3 interpreter sequence and run.

## Convert P-file to 'Localizer.h5'

```
>> pfile = 'P12345.7';
>> readoutFile = 'readout.mod';
>> pfile2hdf(pfile, 1, readoutFile);  
```



