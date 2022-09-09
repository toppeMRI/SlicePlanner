
// HDF5
import ch.systemsx.cisd.hdf5.*;

// Write text file
import java.io.FileWriter;
import java.io.IOException;

public class ROI {

	// ROI id/number 
	int roiId;

	// width, height, and thickness (pixels). TODO: change to cm
	public double w;
	public double h;
	public double t;

	// size limits
	public double wmin;   // minimum dimension (pixels). TODO: change to cm
	public double hmin;
	public double tmin;
	public double wmax;   // pixels. TODO: change to cm
	public double hmax;
	public double tmax;

	// position of center of ROI (pixels), offset from Iso-center
	public double x;
	public double y;
	public double z;

	// 3x3 rotation matrix
	public double[][] rotmat = {   
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	};

	// "Full" class constructor
	ROI(double[] size, double[] center, double[] minSize, double[] maxSize, double[][] rotationMatrix, int id) {

		roiId = id;

		// size
		w = size[0];
		h = size[1];
		t = size[2];

		// center location
		x = center[0];
		y = center[1];
		z = center[2];

		// size limits
		wmin = minSize[0];
		hmin = minSize[1];
		tmin = minSize[2];
		wmax = maxSize[0];
		hmax = maxSize[1];
		tmax = maxSize[2];

		// rotation matrix
		rotmat = rotationMatrix;
	}

	// A convenient "default" constructor
	ROI(int nx, int ny, int nz, int id) {

		roiId = id;

		w = 3.*nx/4.;	
		h = ny/2.;	
		t = nz/8.;	

		x = 0.;
		y = 0.;
		z = 0.;

		wmin = 5.;
		hmin = 5.;
		tmin = 5.;
		wmax = nx;
		hmax = ny;
		tmax = nz;

		rotmat = new double[][]{    
			{1., 0., 0.},
			{0., 1., 0.},
			{0., 0., 1.}
		};
	}

	// Create directly from file
	public ROI(String fname, int roiId) {
		this.roiId = roiId;
		loadFromHDF5(fname);
	}

	// Make it printable
	@Override
	public String toString() {
			String s = "w/h/t=";
			s += Integer.toString((int)Math.round(w)) + "/";
			s += Integer.toString((int)Math.round(h)) + "/";
			s += Integer.toString((int)Math.round(t)) + ", center x/y/z=";
			s += Integer.toString((int)Math.round(-x)) + "/";
			s += Integer.toString((int)Math.round(-y)) + "/";
			s += Integer.toString((int)Math.round(-z));
			s += "\nRotation matrix:\n";
			for(int i=0;i<3;i++) {    
				for(int j=0;j<3;j++) {    
					s += rotmat[i][j] + " ";
				}
				s += "\n";
			}    

			return s;
	}

	// Write this ROI to HDF5 file. If file doesn't already exist it will be created.
	public void appendToHDF5(String fname) {

		IHDF5SimpleWriter writer = HDF5Factory.open(fname);

		int i = roiId + 1;   /* Start indexing at 1 (Matlab-style) */

		writer.writeDouble("ROI" + i + "/dimensions/width", w);
		writer.writeDouble("ROI" + i + "/dimensions/height", h);
		writer.writeDouble("ROI" + i + "/dimensions/thickness", t);

		writer.writeDouble("ROI" + i + "/dimensions/minWidth", wmin);
		writer.writeDouble("ROI" + i + "/dimensions/minHeight", hmin);
		writer.writeDouble("ROI" + i + "/dimensions/minThickness", tmin);

		writer.writeDouble("ROI" + i + "/dimensions/maxWidth", wmax);
		writer.writeDouble("ROI" + i + "/dimensions/maxHeight", hmax);
		writer.writeDouble("ROI" + i + "/dimensions/maxThickness", tmax);

		writer.writeDouble("ROI" + i + "/center/x", -x);
		writer.writeDouble("ROI" + i + "/center/y", -y);
		writer.writeDouble("ROI" + i + "/center/z", -z);

		// write rotation matrix as 1D array, in row-major form (to avoid row/column swaps which HDF5 seems to like to do)
		double[] rot = new double[9];
		for(int ii=0;ii<3;ii++) {    
			for(int jj=0;jj<3;jj++) {    
				rot[ii*3+jj] = rotmat[ii][jj];
			}
		}    
		writer.writeDoubleArray("ROI" + i + "/rotmat", rot);

		// Temporary feature: Also write rotmat.txt which can be loaded directly into toppev3 
		/*
		if (roiId == 0) {
			try {
				FileWriter txtwriter = new FileWriter("rotmat.txt", false);
				for(int l = 0; l<9; l++) {
					txtwriter.write(Double.toString(rot[l]) + "\n");
				}
				txtwriter.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		*/

		// Write (normal) distance from scan plane to isocenter (for calculating RF slice-select frequency offset)
		double[] n = new double[]{0., 0., 1.};   // normal vector along thickness dimension, prior to rotating
		n = Utils.vecrot(rotmat, n);         // normal vector of scan plane
		double[] pIso = new double[]{0., 0., 0};       // iso-center
		double[] pRoiCenter = new double[]{x, y, z};   // center of ROI (offset from iso-center)
		double[] pi = Utils.getRayPlaneIntersection(n, pIso, pRoiCenter, n);  // point of intersection
		double[] piVec = Utils.vecsubtract(pi, pIso);                         // vector from iso-center to intersection point
		double d = Math.sqrt(piVec[0]*piVec[0] + piVec[1]*piVec[1] + piVec[2]*piVec[2]);
		if (Utils.vecdot(n,piVec) > 0) {   /* positive z (in pixel units) = Inferior direction, so need to negate */
			d = -d;
		}
		writer.writeDouble("ROI" + i + "/scanPlaneToIsocenterDistance", d);
		//System.out.println(d);

		System.out.println(this);

		writer.close();
	}

	// Load this ROI from HDF5 file.
	public void loadFromHDF5(String fname) {

		IHDF5SimpleReader reader = HDF5Factory.openForReading(fname);

		int i = roiId + 1;

		w = reader.readDouble("ROI" + i + "/dimensions/width");
		h = reader.readDouble("ROI" + i + "/dimensions/height");
		t = reader.readDouble("ROI" + i + "/dimensions/thickness");

		wmin = reader.readDouble("ROI" + i + "/dimensions/minWidth");
		hmin = reader.readDouble("ROI" + i + "/dimensions/minHeight");
		tmin = reader.readDouble("ROI" + i + "/dimensions/minThickness");

		wmax = reader.readDouble("ROI" + i + "/dimensions/maxWidth");
		hmax = reader.readDouble("ROI" + i + "/dimensions/maxHeight");
		tmax = reader.readDouble("ROI" + i + "/dimensions/maxThickness");

		x = reader.readDouble("ROI" + i + "/center/x");
		y = reader.readDouble("ROI" + i + "/center/y");
		z = reader.readDouble("ROI" + i + "/center/z");

		double[] rot = new double[9];
		rot = reader.readDoubleArray("ROI" + i + "/rotmat");
		for(int ii=0;ii<3;ii++) {    
			for(int jj=0;jj<3;jj++) {    
				rotmat[ii][jj] = rot[ii*3+jj];
			}
		}    

		reader.close();
	}

}
