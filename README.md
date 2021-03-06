# SlicePlanner

GUI for prescribing one (or more) 3D rectangular ROIs off of a stack of 2D images.  

Primary goal is to emulate the slice prescription experience on commercial MRI scanners, 
in a way that is compatible 
with vendor-agnostic MRI pulse programming environments such as 
[Pulseq](http://pulseq.github.io/) and [TOPPE](https://toppemri.github.io).  

Input: The file **Localizer.h5** containing a 3D image volume in HDF5 format. See ./LocalizerScan/.

Output: The file **ROI.h5** containing ROI parameters (size, location, rotation, offset from iso-center).

For Java install+use notes, see ./Doc/.

See also [./GUI/README.md](./GUI/README.md)

![GUI screenshot](Resources/gui.jpg)
