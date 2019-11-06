# Java source code for the SlicePlanner GUI


## Overview

Each 3D ROI is internally represented as a rectangular box, and is encapsulated in the class 'ROI' (**ROI.java**). 
The intersection of this box with an arbitrary 2D viewing plane is calculated on-the-fly. 
The underlying geometric calculations are explicitly “hand-coded” and do not rely on Java libraries, and can be easily reused in other programming (e.g., C++) contexts. 
These elementary calculations are contained in **Utils.java**.

The graphical user interface (GUI) is written in JavaFX. JavaFX-dependent code is confined to **SlicePlanner.java** and **MRImageView.java**.

Using the mouse, the user can interactively resize, move, and rotate each ROI, by manipulating the intersection of the ROI with axial, sagittal, and coronal views of the subject. 

In addition, the center slice within each ROI is displayed in separate windows, as an additional guide to ROI placement (a feature not commonly found on commercial scanners). The GUI also supports slice scrolling with the mouse wheel. 

Clicking the ‘Export ROIs’ button writes all ROIs to a single HDF5 file that can be loaded into Matlab with the ‘toppe.getroi()’ function in the TOPPE [1] Matlab toolbox (https://toppemri.github.io) for subsequent use in the design of vendor-agnostic sequences based on TOPPE or Pulseq. 

This graphical tool can be run on commonly used personal computers, and does not rely on graphics processing units (GPU) or other advanced graphics features.


## Starting the GUI (Linux)

Prerequisite: The file **Localizer.h5** must exist in the local path. 

To start the GUI:
```
$ ./start
```

