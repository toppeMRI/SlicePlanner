//package Utils;

import java.util.Arrays;

public class Utils {
	
	// Get intersection point between a line segment and a plane (in 3D space)
	// See https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
	//
	// Inputs:
   //   l0   start point of line segment
   //   l1   end point of line segment
	//   p0   a point on the 2D plane
	//   n    normal vector of the plane
	//
	//	Output:
	//   intersection point, double[3]. Returns null if no intersection point found.
	//
	public static double[] getLinePlaneIntersection(double[] l0, double[] l1, double[] p0, double[] n) {

		double[] l = vecsubtract(l1,l0);    // A vector pointing along the direction of the line.

		double[] pi = getRayPlaneIntersection(l, l0, p0, n);   // intersection point

		// does pi lie on line segment l?
		if (pi != null) {
			double minx = minarr(new double[]{l0[0],l1[0]});
			double maxx = maxarr(new double[]{l0[0],l1[0]});
			double miny = minarr(new double[]{l0[1],l1[1]});
			double maxy = maxarr(new double[]{l0[1],l1[1]});
			double minz = minarr(new double[]{l0[2],l1[2]});
			double maxz = maxarr(new double[]{l0[2],l1[2]});
			if (pi[0]>=minx && pi[0]<=maxx && pi[1]>=miny && pi[1]<=maxy && pi[2]>=minz && pi[2]<=maxz) {
				// keep it
			} else {
				// Intersection point does not lie on the line (segment)
				pi = null;
			}
		}

		return pi;
	}

		
	// Get intersection point between a ray (of infinite length) and a 2D plane (in 3D space)
	// See https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
	// Inputs:
	//  l     vector pointing along ray 
  	//  l0    a point on the ray
	//  p0    a point on the 2D plane
	//  n     normal vector of the plane
	//
	public static double[] getRayPlaneIntersection(double[] l, double[] l0, double[] p0, double[] n) {

		double[] pi = new double[3];        // Intersection point

		double tol = 1e-16;

		if (Math.abs(vecdot(l, n)) < tol) {
			if (Math.abs(vecdot(vecsubtract(p0,l0),n)) == 0.0) {
				// l is in the viewing plane, so we'll include both endpoints of l in pi
				// actually, just ignore for now -- fix later?
				//pi[cnt++] = p[pairs[i][0]];
				//pi[cnt++] = p[pairs[i][1]];
				//System.out.println("l in viewing plane");
				pi = null;  // ok in in the context of this GUI
			}
			else {
				// No intersection found: l is parallel with, but not in, the viewing plane
				pi = null;
			}
			pi = null;
		} 
		else {
			// There is a point of intersection between an infinite ray along the line and the 2D plane.
			double d = vecdot(vecsubtract(p0,l0),n)/vecdot(l,n);
			for (int k = 0; k < 3; k++) {
				pi[k] = d*l[k] + l0[k];
			}
		}

		return pi;
	}

	// Trilinear interpolation
	// Uses notation from https://en.wikipedia.org/wiki/Trilinear_interpolation
	public static double trilinearInterp(double x0, double x1, double y0, double y1, double z0, double z1, 
		double c000, double c100, double c010, double c110, double c001, double c101, double c011, double c111, 
		double x, double y, double z) { 

		double c0 = bilinearInterp(x0, x1, y0, y1, c000, c100, c010, c110, x, y);
		double c1 = bilinearInterp(x0, x1, y0, y1, c001, c101, c011, c111, x, y);

		double zd = (z-z0)/(z1-z0);

		double val = c0*(1-zd) + c1*zd;

		return val;
	}

	public static double bilinearInterp(double x0, double x1, double y0, double y1, 
		double c00, double c10, double c01, double c11, 
		double x, double y) {

		double c0 = linearInterp(x0, x1, c00, c10, x);
		double c1 = linearInterp(x0, x1, c01, c11, x);

		double yd = (y-y0)/(y1-y0);

		double val = c0*(1-yd) + c1*yd;
		
		return val;
	}

	public static double linearInterp(double x0, double x1, double c0, double c1, double x) {
		double xd = (x-x0)/(x1-x0);
		double val = c0*(1-xd) + c1*xd;
		return val;
	}


	// Get 2D 'slice' of a 3D image volume. 
	// Viewing plane is determined by the given rotation matrix and translation (offset).
	//
   // Inputs:
	//   ims      int[nx][ny][nz]    image volume
   //   rotmat   double[3][3]       rotation matrix
	//   offset   double[3]          offset of center of viewing plane to iso-center (pixels)
	//
   // Output:
   //   imObl    int[nx][ny]        2D slice through ims   
	//
	public static int[][] getSlice(int[][][] ims, double[][] rotmat, double[] offset) {

		int nx = ims.length;
		int ny = ims[0].length;
		int nz = ims[0][0].length;
	
		// Construct grid 
		int xmin = -nx/2;
		int ymin = -ny/2;
		int zmin = -nz/2;
		int xmax = nx/2-1;
		int ymax = ny/2-1;
		int zmax = nz/2-1;

		double[][] X = new double[nx][ny];
		double[][] Y = new double[nx][ny];
		double[][] Z = new double[nx][ny];

		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				X[i][j] = i;
				Y[i][j] = j;
				Z[i][j] = nz/2.-0.5; // by definition, the initial (unrotated, untranslated) viewing plane lies in the xy plane
			}
		}

		// Rotate and translate grid points
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				double[] p = new double[]{X[i][j]-nx/2.-0.5, Y[i][j]-ny/2.-0.5, 0.};
				p = vecrot(rotmat, p);
				X[i][j] = nx/2. + 0.5 + p[0] + offset[0];
				Y[i][j] = ny/2. + 0.5 + p[1] + offset[1];
				Z[i][j] += p[2] + offset[2];
			}
		}

		// Interpolate onto viewing plane grid
		int[][] imObl = new int[nx][ny];
		int x0, y0, z0, x1, y1, z1;
		double c000, c100, c010, c110, c001, c101, c011, c111;
		double val;
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				x0 = (int) Math.floor(X[i][j]);
				y0 = (int) Math.floor(Y[i][j]);
				z0 = (int) Math.floor(Z[i][j]);
				if (x0>0 && y0>0 && z0>0 && x0<nx-1 && y0<ny-1 && z0<nz-1) {
					x1 = x0 + 1;
					y1 = y0 + 1;
					z1 = z0 + 1;
					c000 = (double) ims[x0][y0][z0];
					c100 = (double) ims[x1][y0][z0];
					c010 = (double) ims[x0][y1][z0];
					c110 = (double) ims[x1][y1][z0];
					c001 = (double) ims[x0][y0][z1];
					c101 = (double) ims[x1][y0][z1];
					c011 = (double) ims[x0][y1][z1];
					c111 = (double) ims[x1][y1][z1];
					val = Utils.trilinearInterp(x0, x1, y0, y1, z0, z1, 
						c000, c100, c010, c110, c001, c101, c011, c111, 
						X[i][j], Y[i][j], Z[i][j]);
					imObl[i][j] = (int) Math.round(val);
				}
				else {
					// grid point lies outside the image volume
					imObl[i][j] = (int) 0;
				}
			}
		}

		return imObl;
	}
		
	// Convert 2D grayscale (int) image to 1D RGB (byte) array suitable for passing to PixelWriter.setPixels
	// The first dimension is the "fast" dimension and corresponds to "x" in the Java scene (i.e., left/right)
	public static byte[] grayscaleIm2rgb(int[][] im, int nx, int ny) {

		byte[] imRgb = new byte[3*nx*ny]; 

		int cnt = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				imRgb[cnt]   = (byte) im[i][j];
				imRgb[cnt+1] = (byte) im[i][j];
				imRgb[cnt+2] = (byte) im[i][j];
				cnt += 3;
			}
		}

		return imRgb;

	}

	// Return the (clockwise positive) angle from 0 - 360 degrees, 
	// between two lines defined by points (x,y) and (px,py)
	public static double clockAngle (double x, double y, double px, double py) {
		double dx = x - px;
		double dy = y - py;

		double angle = Math.abs(Math.toDegrees(Math.atan2(dy, dx)));

		if(dy < 0) {
			angle = 360 - angle;
		}

		//angle = 360 - angle;   // return counter-clockwise angle

		return angle;
	}

	// Multiply two rotation matrices
	public static double[][] rotmult(double[][] a, double[][] b) {
		double c[][] = new double[3][3]; 
   
		for(int i=0;i<3;i++) {    
			for(int j=0;j<3;j++) {    
				c[i][j] = 0;      
				for(int k=0;k<3;k++) {
					c[i][j] += a[i][k]*b[k][j];      
				}
			}
		}    

		return c;
	} 

	// Apply rotation matrix 'rot' to a point 'b' in 3D space
	public static double[] vecrot(double[][] rot, double[] b) {
		double c[] = new double[3]; 
   
		for(int i=0;i<3;i++) {    
			c[i]=0;      
			for(int k=0;k<3;k++) {
				c[i] += rot[i][k]*b[k];      
			}
		}    

		return c;
	} 

	// Dot product of two points (vectors) in 3D space
	public static double vecdot(double[] a, double[] b) {
		double c = 0;
   
		for(int k=0;k<3;k++) {
			c += a[k]*b[k];      
		}

		return c;
	} 

	// Subtract two points (vectors) in 3D space
	public static double[] vecsubtract(double[] a, double[] b) {
		double[] c = new double[3];
   
		for(int k=0;k<3;k++) {
			c[k] = a[k]-b[k];      
		}

		return c;
	} 

	// Get max value in a 1D array
	public static double maxarr(double[] arr) {
		double[] a = new double[arr.length];
		for (int i=0; i < arr.length; i++) {
			a[i] = arr[i];
		}
		Arrays.sort(a);
		return a[a.length-1];
	}

	// Get min value in a 1D array
	public static double minarr(double[] arr) {
		double[] a = new double[arr.length];
		for (int i=0; i < arr.length; i++) {
			a[i] = arr[i];
		}
		Arrays.sort(a);
		return a[0];
	}

	public static void printPoints(double[][] pi, int cnt) {
		String s = ""; //"pi[1:cnt][:] = ";
		for (int i = 0; i < cnt; i++) {
			s += "(";
			for (int j = 0; j < 3; j++) {
				s += Integer.toString((int) Math.round(pi[i][j]));
				if (j < 2) {
					s += ",";
				}
			}
			s += ")";
		}
		System.out.println(s);
	}

	public static void printPoint(double[] pi) {
		String s = "("; //"pi[1:cnt][:] = ";
		for (int j = 0; j < 3; j++) {
			s += Integer.toString((int) Math.round(pi[j]));
			if (j < 2) {
				s += ",";
			}
		}
		s += ")";
		System.out.println(s);
	}

	public static void printIntVec(int[] p, int cnt) {
		String s = "";
		s += "(";
		for (int i = 0; i < cnt-1; i++) {
			s += p[i] + ",";
		}
		s += p[cnt-1] + ")";
		System.out.println(s);
	}

	public static void printDoubleVec(double[] p, int cnt) {
		String s = "";
		s += "(";
		for (int i = 0; i < cnt-1; i++) {
			s += p[i] + ",";
		}
		s += p[cnt-1] + ")";
		System.out.println(s);
	}

	// Keep a few rotation matrices handy (for defining fixed sag/axi/cor views using Utils.getSlice())
	public static double rotmatAxi[][] = {   
		{1., 0., 0.},
		{0., 1., 0.},
		{0., 0., 1.}
	};
	public static double rotmatPlusY[][] = {   
		{ 0., 0., 1.},
		{ 0., 1., 0.},
		{-1., 0., 0.}
	};
	public static double rotmatMinusX[][] = {   
		{1., 0.,  0.},
		{0., 0., -1.},
		{0., 1.,  0.}
	};
	public static double rotmatSag[][] = rotmult(rotmatMinusX, rotmatPlusY);
}
