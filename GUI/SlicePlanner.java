/*
Graphical slice prescription off of 3-plane (axial, sagittal, coronal) views of a 3D image volume.

Input:
   3D MRI image volume (stack of 2D slices)
	
Output:
   3D rectangular ROI(s) (see ROI.java) 

Developed in openjfx-12.0.2 and openjdk-12.0.2
*/

import javafx.stage.Stage; 

import javafx.scene.Node;

import javafx.application.Application; 
import javafx.scene.Scene; 
import javafx.scene.paint.Color; 

import ch.systemsx.cisd.hdf5.*;

import javafx.scene.transform.Transform;
import javafx.scene.transform.Rotate; 
import javafx.scene.transform.Scale; 
import javafx.scene.transform.Translate; 

import javafx.scene.shape.Line;
import javafx.scene.shape.Circle;
import javafx.scene.shape.StrokeType;
import javafx.scene.shape.Rectangle; 

import javafx.scene.layout.GridPane;
import javafx.scene.layout.Region;

import javafx.scene.input.InputEvent;
import javafx.scene.input.MouseEvent;
import javafx.event.EventHandler;
import javafx.scene.Cursor;

import javafx.geometry.Insets;

import javafx.scene.image.ImageView;  
import javafx.scene.image.WritableImage;
import javafx.scene.image.PixelWriter;
import javafx.scene.image.PixelFormat;

import java.io.FileInputStream;
import java.io.File;

import java.nio.ByteBuffer;

import javafx.animation.AnimationTimer;

import java.util.Arrays;

import javafx.scene.control.TextField;
import javafx.scene.control.Label;
import javafx.scene.layout.VBox;

import javafx.event.ActionEvent;

import javafx.scene.input.ScrollEvent;
import javafx.geometry.Bounds;

import javafx.scene.control.Button;

import javafx.scene.shape.Polyline;
import javafx.scene.shape.Polygon;

import java.io.FileWriter;
import java.io.IOException;

import javafx.scene.control.CheckBox;

public class SlicePlanner extends Application {

	static int nROI = 3;              // max number of ROIs to draw in this application.

	static ROI[] rois = new ROI[nROI];   // list of ROIs

	static int nx, ny, nz;               // image matrix size
	static double isox, isoy, isoz;      // iso-center of image volume (in matrix units, i.e., center = (nx/2,ny/2,nz/2))

	static double zoom = 1.3;            // overall app zoom 

	static int[][][] ims;                // 3D image volume (scaled o [0,255])

	static CheckBox[] cb = new CheckBox[nROI];

	MRImageView imageViewAxi;
	MRImageView imageViewSag;
	MRImageView imageViewCor;
	MRImageView[] imageViewObl;

	SelectionBox[] sbAxi;
	SelectionBox[] sbSag;
	SelectionBox[] sbCor;

	GridPane grid;
	GridPane controlGrid;

	// Load ROIs from file 
	public static class LoadRoiButton extends Button {

		public LoadRoiButton(String label) {
			this.setText(label);

			setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent e) {
					loadFile("ROI.h5");   // fill rois[] array with info from file
					//updateGUI();          // update checkboxes
				}
			});
		}

		private void loadFile(String fname) {
			//IHDF5SimpleReader reader = HDF5Factory.openForReading(fname);
			//nROI = reader.readInt("nROI");
			//reader.close();

			for (int i = 0; i < nROI; i++) {
				rois[i] = new ROI(fname, i);
			}
		}
	}

	// Save ROIs to file
	public static class SaveRoiButton extends Button {

		public SaveRoiButton(String label) {

			this.setText(label);

			setOnAction(new EventHandler<ActionEvent>() {
				@Override
				public void handle(ActionEvent e) {
					writeFile("ROI.h5");
				}
			});
		}

		public void writeFile(String fname) {

			// Delete HDF5 file if it already exists
			File file = new File(fname);
			file.delete();
			
			// Append roi to file
			for (int i = 0; i < nROI; i++) {
				rois[i].appendToHDF5(fname);
			}

			//IHDF5SimpleWriter writer = HDF5Factory.open(fname);
			//writer.writeInt("nROI", nROI);
			//writer.close();
		}
	}

	// Main entry point of application
	@Override
	public void start(Stage stage) throws Exception {

		// Load image volume from HDF5 file 
		IHDF5SimpleReader reader = HDF5Factory.openForReading("Localizer.h5");
		
		nx = reader.readInt("/Dims/nx");
		ny = reader.readInt("/Dims/ny");
		nz = reader.readInt("/Dims/nz");
		isox = nx/2.;
		isoy = ny/2.;
		isoz = nz/2.;

		ims = new int[ny][nx][nz];

		for (int k = 0; k < nz; k++) {
			int[][] imAxi = reader.readIntMatrix("/Ax/slice" + Integer.toString(k+1));
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					// NB! IHDF5SimpleReader seems to return the transpose of the 2D matrix in the hdf5 file.
					// For displaying slices, scene x = first matrix dimension (see Utils.grayscaleIm2rgb)
					ims[j][i][k] = imAxi[i][j];  // ims is now consistent with the matrix written to hdf5 file in Matlab (scene x = first dimension of Matlab matrix 'imsos')
				}
			}
		}

		reader.close();

		// Initialize the ROI(s).   TODO: read initial ROIs from config file
		for (int id = 0; id < nROI; id++) {		
			double[] size = new double[]{3.*nx/4, 3.*nx/4, nz/8.};
			double[] center = new double[]{0., 0., 0.}; // -nz/4 + 40*id};      // offset from iso-center
			double[] minSize = new double[]{20., 5., 5.};
			double[] maxSize = new double[]{nx, ny, nz};
			rois[id] = new ROI(size, center, minSize, maxSize, Utils.rotmatAxi, id);
			//rois[id] = new ROI(nx, ny, nz, id);
		}

		// Create MRImageView objects for displaying 2D slices
		imageViewAxi = new MRImageView(nx, ny, -nz/2.+1, nz/2.-1);
		imageViewSag = new MRImageView(ny, nz, -nx/2.+1, nx/2.-1);
		imageViewCor = new MRImageView(nz, nz, -ny/2.+1, ny/2.-1);
		imageViewObl = new MRImageView[nROI];
		for (int i = 0; i < nROI; i++) {
			imageViewObl[i] = new MRImageView(ny,nx);
		}

		//
		// GridPane
		//
		grid = new GridPane();
		grid.setPadding(new Insets(10, 10, 10, 10));
		int viewsize = (int)(zoom*Math.max(nx, ny));
		grid.setMinSize(3*viewsize + 400, (int)Math.round(2*viewsize + 200));
		grid.setVgap(30);
		grid.setHgap(30);
		grid.getTransforms().add(new Scale(zoom, zoom, 0, 0));    // Zoom in on entire scene

		// Control panel
		controlGrid = new GridPane();
		controlGrid.setPadding(new Insets(10, 10, 10, 10));
		controlGrid.setVgap(30);
		for (int i = 0; i < nROI; i++) {
			cb[i] = new CheckBox("ROI " + (i+1));
			cb[i].setSelected(true);
			controlGrid.add(cb[i], 0, i);
		}

		LoadRoiButton loadRoiButton = new LoadRoiButton("Load ROIs");
		controlGrid.add(loadRoiButton, 0, nROI);

		SaveRoiButton saveRoiButton = new SaveRoiButton("Export ROIs");
		controlGrid.add(saveRoiButton, 0, nROI+1);

		grid.add(controlGrid, 0, 0);

		// Axial view
		grid.add(imageViewAxi, 1, 0);
		sbAxi = new SelectionBox[nROI];
		for (int i = nROI-1; i >= 0; i--) {   // add in reverse order so ROI 1 is on top
			sbAxi[i] = new SelectionBox("Axi", new double[]{0., 0., 1.}, i);
			grid.add(sbAxi[i], 1, 0);
		}

		// Sagittal view
		grid.add(imageViewSag, 2, 0);
		sbSag = new SelectionBox[nROI];
		for (int i = nROI-1; i >= 0; i--) {
			sbSag[i] = new SelectionBox("Sag", new double[]{1., 0., 0.}, i);
			grid.add(sbSag[i], 2, 0);
		}

		// Coronal view
		grid.add(imageViewCor, 3, 0);
		sbCor = new SelectionBox[nROI];
		for (int i = nROI-1; i >= 0; i--) {
			sbCor[i] = new SelectionBox("Cor", new double[]{0., 1., 0.}, i);
			grid.add(sbCor[i], 3, 0);
		}

		// Oblique view(s). Shows the logical z=0 slice (width/height plane) of each ROI.
		for (int i = 0; i < nROI; i++) {
			grid.add(imageViewObl[i], i+1, 1);
		}

		// Create a new scene
		Scene scene = new Scene(grid);
		stage.setScene(scene);
		stage.setTitle("SlicePlanner");
		//stage.setMinHeight(600);
		//stage.setMinWidth(1200);
		stage.show();

		// Create animation timer so view refreshes periodically. This is a really nice feature of JavaFX.
		int[][] imzero = new int[nx][ny];
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				imzero[i][j] = 200;
			}
		}

		new AnimationTimer() {
			private long lastUpdate;
			@Override
			public void handle(long now) {
				if (now - lastUpdate > 10_000_000) {
					lastUpdate = now;

					int[][] im;
					double[] offset;

					// Update viewports
					
					// Axial view
					offset = new double[]{0., 0., imageViewAxi.sliceOffset};
					im = Utils.getSlice(ims, Utils.rotmatAxi, offset);
					imageViewAxi.setPixels(im);
					for (int i = 0; i < nROI; i++) {
						sbAxi[i].offset = imageViewAxi.sliceOffset;
						sbAxi[i].draw();
					}

					// Sagittal view
					offset = new double[]{imageViewSag.sliceOffset, 0., 0.};
					im = Utils.getSlice(ims, Utils.rotmatSag, offset);
					imageViewSag.setPixels(im);
					for (int i = 0; i < nROI; i++) {
						sbSag[i].offset = imageViewSag.sliceOffset;
						sbSag[i].draw();
					}

					// Coronal view
					offset = new double[]{0., imageViewCor.sliceOffset, 0.};
					im = Utils.getSlice(ims, Utils.rotmatMinusX, offset);
					imageViewCor.setPixels(im);
					for (int i = 0; i < nROI; i++) {
						sbCor[i].offset = imageViewCor.sliceOffset;
						sbCor[i].draw();
					}

					// Oblique views through the center of each ROI
					// First, paint all gray.
					for (int i = 0; i < nROI; i++) {
						if (cb[i].isSelected()) {
							offset = new double[]{rois[i].x, rois[i].y, rois[i].z};
							im = Utils.getSlice(ims, rois[i].rotmat, offset);
							imageViewObl[i].setPixels(im);
						} 
						else {
							imageViewObl[i].setPixels(imzero);
						}
					}
				}
			}
		}.start();
	}

	public static void updateGUI() {

		// Set ROI checkboxes
		for (int i = nROI; i < nROI; i++) {
			cb[i].setSelected(false);
			//controlGrid.getChildren().remove(cb[i]);
		}
		for (int i = 0; i < nROI; i++) {
			cb[i].setSelected(true);
		}
	}

	// 2D polygon outline of ROI in viewing plane.
	// Can be rotated, dragged, and reshaped -- the underlying ROI is updated accordingly.
	public static class SelectionBox extends Region {
		// Based on https://stackoverflow.com/questions/12644031/drag-mouse-to-rotate-and-scale-in-java

		String orientation;  // "Axi", "Sag", or "Cor"

		double offset = 0.;                   // offset of this view plane from iso-center (pixels)

		int roiNum = 0;                      // index into arrays of ROIs

		// Defines viewing orientation. 
		double[] normalVec = new double[3];   

		double angle;                         // in-plane rotation angle

		Polygon polygon = new Polygon();

		private enum Position {
			Top, BottomRight, Right, Bottom, BottomLeft, Left; 
		}

		// Create the corners
		private Rectangle rb;

		// Size of corner boxes
		private double cornerSize = 8;

		// Create a new rotate transform
		private final Rotate rotate = new Rotate();
		{
			//getTransforms().add(rotate);
			rotate.setPivotX(cornerSize);
			rotate.setPivotY(cornerSize);
		}

		// Circle which is dragged to rotate the box
		private final Circle rotateCircle;

		// Box at center for moving the SelectionBox
		private Rectangle tb;    // center box

		// Variables to store mouse x and y
		private double x, y;

		public SelectionBox (String orient, double[] nv, int roiNum) {

			orientation = orient;

			for (int i = 0; i < 3; i++) {
				normalVec[i] = nv[i];
			}

			this.roiNum = roiNum;

			setPickOnBounds(false);    // needed to make MRImageView node receive mouse (scroll) events

			// Create the circle which can be dragged to rotate the box
			rotateCircle = new Circle(5);
			rotateCircle.setFill(Color.PINK);
			rotateCircle.setStroke(Color.rgb(0,0,0, 0.75));

			// Make it draggable
			rotateCircle.addEventHandler(MouseEvent.MOUSE_PRESSED, new EventHandler<MouseEvent>() {
				@Override public void handle(MouseEvent event) {
					setMouse(event.getSceneX(), event.getSceneY());
				}
			});

			// When it's dragged, rotate the box and update the ROI rotation matrix
			rotateCircle.addEventHandler(MouseEvent.MOUSE_DRAGGED, new EventHandler<MouseEvent>() {
				@Override public void handle(MouseEvent event) {

					// Used to get the scene position of the corner of the box
					Transform localToScene = getLocalToSceneTransform();

					double x1 = getMouseX();
					double y1 = getMouseY();

					double x2 = event.getSceneX();
					double y2 = event.getSceneY();

					double px = rotate.getPivotX() + localToScene.getTx();
					double py = rotate.getPivotY() + localToScene.getTy();

					// Work out the angle rotated 
					double th1 = Utils.clockAngle(x1, y1, px, py);
					double th2 = Utils.clockAngle(x2, y2, px, py);

					// keep track of net rotation angle in this view
					angle += th2 - th1;
					if (angle > 360) {
						angle -= 360;
					}
					if (angle < 0) {
						angle += 360;
					}
					//rotate.setAngle(angle);

					// Update mouse coordinates
					setMouse(event.getSceneX(), event.getSceneY());

					// Update ROI rotation matrix
					double phi = -Math.toRadians(th2 - th1);   // Note negation
					double cphi = Math.cos(phi);
					double sphi = Math.sin(phi);

					if (orientation.equals("Axi")) {
						double rotmat[][] = {   
							{ cphi, sphi, 0.},           // sign of sphi terms chosen to be consisten with EPIC
							{-sphi, cphi, 0.},
							{   0.,   0., 1.}
						};
						rois[roiNum].rotmat = Utils.rotmult(rotmat, rois[roiNum].rotmat);
					}

					if (orientation.equals("Sag")) {
						double rotmat[][] = {   
							{1.,   0.,    0.},
							{0., cphi, sphi}, 
							{0., -sphi,  cphi}
						};
						rois[roiNum].rotmat = Utils.rotmult(rotmat, rois[roiNum].rotmat);
					}

					if (orientation.equals("Cor")) {
						double rotmat[][] = {   
							{ cphi, 0., sphi},
							{   0., 1.,   0.},
							{-sphi, 0., cphi}
						};
						rois[roiNum].rotmat = Utils.rotmult(rotmat, rois[roiNum].rotmat);
					}
				
					draw();
				}
			});

			// Build the corners
			rb = buildResizeBox(0,0, Position.BottomRight);

			// Build the translation box
			tb = buildTranslationBox(0,0);

			// Display cursor when mouse hovers over control boxes
			rb.setCursor(Cursor.CROSSHAIR);
			tb.setCursor(Cursor.CROSSHAIR);
			rotateCircle.setCursor(Cursor.CROSSHAIR);

			getChildren().addAll(polygon, rotateCircle, tb, rb);

		} // end of SelectionBox constructor

		// Draw the selection box
		public void draw () {
			// Figure out intersection between ROI (box) edges and viewing plane

			double[][] p = new double[8][3];    // corner points

			// Corner points (centered at 0, before rotating)
			p[0][0] = -rois[roiNum].w/2;
			p[0][1] = -rois[roiNum].h/2;
			p[0][2] = -rois[roiNum].t/2;
			p[1][0] = -rois[roiNum].w/2;
			p[1][1] =  rois[roiNum].h/2;
			p[1][2] = -rois[roiNum].t/2;
			p[2][0] =  rois[roiNum].w/2;
			p[2][1] =  rois[roiNum].h/2;
			p[2][2] = -rois[roiNum].t/2;
			p[3][0] =  rois[roiNum].w/2;
			p[3][1] = -rois[roiNum].h/2;
			p[3][2] = -rois[roiNum].t/2;

			p[4][0] = -rois[roiNum].w/2;
			p[4][1] = -rois[roiNum].h/2;
			p[4][2] =  rois[roiNum].t/2;
			p[5][0] = -rois[roiNum].w/2;
			p[5][1] =  rois[roiNum].h/2;
			p[5][2] =  rois[roiNum].t/2;
			p[6][0] =  rois[roiNum].w/2;
			p[6][1] =  rois[roiNum].h/2;
			p[6][2] =  rois[roiNum].t/2;
			p[7][0] =  rois[roiNum].w/2;
			p[7][1] = -rois[roiNum].h/2;
			p[7][2] =  rois[roiNum].t/2;

			// Rotate and translate corner points
			for (int i = 0; i < 8; i++) {
				p[i] = Utils.vecrot(rois[roiNum].rotmat, p[i]);
				p[i][0] += rois[roiNum].x + isox;
				p[i][1] += rois[roiNum].y + isoy;
				p[i][2] += rois[roiNum].z + isoz;
			}

			// Vectors l representing each edge, and points l0 on each edge
			double[][]  l = new double[12][3];
			double[][] l0 = new double[12][3];

			int[][] pairs = {
				{0, 1},{1, 2},{2, 3},{3, 0},    // front face
				{0, 4},{1, 5},{2, 6},{3, 7},    
				{4, 5},{5, 6},{6, 7},{7, 4}
			};
			
			for (int i = 0; i < 12; i++) {
				l[i] = Utils.vecsubtract(p[pairs[i][0]], p[pairs[i][1]]);
				l0[i] = p[pairs[i][0]];
			}

			// A point p0 in the viewing plane, and normal vector n of the viewing plane
			double[] n = normalVec; 
			double[] p0 = new double[3];
			p0[0] = (isox+offset)*n[0];
			p0[1] = (isoy+offset)*n[1];
			p0[2] = (isoz+offset)*n[2];

			// Intersection points pi between ROI edges and viewing plane (axial/sagittal/coronal)
			double[][] pi = new double[6][3];    // Anywhere from 3 to 6 intersection points between a line and a 3D box

			int cnt = 0;   // counts the number of intersections we've found along the way (between edge and viewing plane)
			for (int i = 0; i < 12; i++) {
				double[] pc = Utils.getLinePlaneIntersection(p[pairs[i][0]], p[pairs[i][1]], p0, n);
				if (pc != null) {
					pi[cnt++] = pc;
				}
			}

			if (cnt < 3 || !cb[roiNum].isSelected()) { 
				// Box does not intersect viewing plane (otherwise there are either 3, 4, 5, or 6 intersection points)
				polygon.setStroke(null);
				polygon.setFill(null);
				rb.setStroke(null);
				rb.setFill(null);
				tb.setStroke(null);
				tb.setFill(null);
				rotateCircle.setStroke(null);
				rotateCircle.setFill(null);
			}
			else {
				//polygon.setStroke(Color.RED);
				polygon.setStrokeWidth(2.0);
				polygon.setFill(null);
				rb.setFill(Color.RED);
				tb.setFill(Color.BLUE);
				rotateCircle.setFill(Color.PINK);

				// Project intersection points onto viewing plane
				double[] xv = new double[cnt];    // here 'xv' and 'yv' are defined to lie in the viewing plane
				double[] yv = new double[cnt];

				for (int i = 0; i < cnt; i++) {
					if (orientation.equals("Axi")) {
						xv[i] = pi[i][0];   // x coordinate
						yv[i] = pi[i][1];   // y coordinate
					}
					else if (orientation.equals("Sag")) {
						xv[i] = pi[i][1];  // y coordinate
						yv[i] = pi[i][2];  // z coordinate
					}
					else {               // "Cor"
						xv[i] = pi[i][0];  // x coordinate
						yv[i] = pi[i][2];  // z coordinate
					}
				}
		
				// Before connecting the points with lines, need to sort them in clockwise order
				// We will do this by calculating angle of radial line from center of intersection points and sorting based on that.
				double xvc = (Utils.minarr(xv)+Utils.maxarr(xv))/2;  // center of the intersection points (in viewing plane)
				double yvc = (Utils.minarr(yv)+Utils.maxarr(yv))/2;

				double[] angles = new double[cnt];
				double[] anglesSorted = new double[cnt];
				for (int i = 0; i < cnt; i++) {
					angles[i] = Utils.clockAngle(xv[i], yv[i], xvc, yvc) ;
					anglesSorted[i] = angles[i];
				}
				Arrays.sort(anglesSorted);

				// Sort xv/yv and put result in xvs/yvs
				double[] xvs = new double[cnt];
				double[] yvs = new double[cnt];
				Double[] pvs = new Double[2*cnt];   // for setting polyline points
				int j2 = 0;
				for (int i = 0; i < cnt; i++) {
					for (int j = 0; j < cnt; j++) {
						if (angles[j] == anglesSorted[i]) {
							xvs[i] = xv[j];
							yvs[i] = yv[j];
							pvs[j2] = new Double(xvs[i]);
							pvs[j2+1] = new Double(yvs[i]);
							j2 += 2;
						}
					}
				}

				// Draw polygon
				polygon.getPoints().setAll(pvs);
				polygon.setStroke(Color.RED);

				// TODO: draw line/arrow showing through-slice (thickness) dimension

				// Place resize box
				// Choose corner placement depending on how much in-plane rotation has occurred in this view
				// TODO: this doesn't work quite right yet
				double a = angle + 45;
				if (a > 360) {
					a -= 360;
				}
				if (a < 0) {
					a += 360;
				}
				if (a >= 0 && a < 90) {
					rb.setX(Utils.maxarr(xv));
					rb.setY(Utils.maxarr(yv));
				} else if (a >= 90 && a < 180) {
					rb.setX(Utils.minarr(xv));
					rb.setY(Utils.maxarr(yv));
				} else if (a >= 180 && a < 270) {
					rb.setX(Utils.minarr(xv));
					rb.setY(Utils.minarr(yv));
				} else {
					rb.setX(Utils.maxarr(xv));
					rb.setY(Utils.minarr(yv));
				}
				//rb.setX(getMouseX()/zoom);
				//rb.setY(getMouseY()/zoom);

				// Place translation box
				tb.setX(Utils.minarr(xv));
				tb.setY(Utils.minarr(yv));

				// Place rotateCircle
				rotateCircle.setTranslateX(Utils.maxarr(xv)+20); // + rotateCircle.getRadius());
				rotateCircle.setTranslateY(Utils.minarr(yv));

				// Set rotate pivot to center of SelectionBox
				rotate.setPivotX((Utils.minarr(xv)+Utils.maxarr(xv))/2);
				rotate.setPivotY((Utils.minarr(yv)+Utils.maxarr(yv))/2);
			}

		}    // End of draw()

		// Save mouse coordinates
		private void setMouse(double x, double y) {
			this.x = x;
			this.y = y;
		}

		private double getMouseX () {
			return x;
		}

		private double getMouseY () {
			return y;
		}

		// Build the resize box (rb) 
		private Rectangle buildResizeBox (double x, double y, final Position pos) {

			// Create the rectangle
			Rectangle r = new Rectangle();
			r.setX(x);
			r.setY(y);
			r.setWidth(cornerSize);
			r.setHeight(cornerSize);

			// Make it draggable
			r.addEventHandler(MouseEvent.MOUSE_PRESSED, new EventHandler<MouseEvent>() {
				@Override public void handle(MouseEvent event) {
					setMouse(event.getSceneX(), event.getSceneY()); // Saves (x,y) coordinates of where mouse is clicked in (this.x, this.y)
				}
			});

			r.addEventHandler(MouseEvent.MOUSE_DRAGGED, new EventHandler<MouseEvent>() {
				@Override public void handle(MouseEvent event) {

					// Get the mouse deltas. Account for overall scene scaling (zoom).
					double dx = (event.getSceneX() - getMouseX())/zoom;  // getMouseX() returns this.x
					double dy = (event.getSceneY() - getMouseY())/zoom;

					setMouse(event.getSceneX(), event.getSceneY());

					// Unit vectors along the width, height, and thickness dimensions
					double[] wv = new double[]{1., 0., 0.};
					wv = Utils.vecrot(rois[roiNum].rotmat, wv);
					double[] hv = new double[]{0., 1., 0.};
					hv = Utils.vecrot(rois[roiNum].rotmat, hv);
					double[] tv = new double[]{0., 0., 1.};
					tv = Utils.vecrot(rois[roiNum].rotmat, tv);

					// Update ROI dimensions
					double[] sv = new double[3]; // Vector representing the stretch we're applying in this view

					if (orientation.equals("Axi")) {
						sv = new double[]{dx,dy,0.};
					}
					else if (orientation.equals("Sag")) {
						sv = new double[]{0.,dx,dy};
					}
					else {
						sv = new double[]{dx,0.,dy};
					}
					double dw = 2 * Utils.vecdot(sv, wv);    
					double dh = 2 * Utils.vecdot(sv, hv);    
					double dt = 2 * Utils.vecdot(sv, tv);    
				
					rois[roiNum].w += dw;
					rois[roiNum].h += dh;
					rois[roiNum].t += dt;

					// Don't let dimensions exceed limits
					rois[roiNum].w = Math.max(rois[roiNum].w, rois[roiNum].wmin);
					rois[roiNum].h = Math.max(rois[roiNum].h, rois[roiNum].hmin);
					rois[roiNum].t = Math.max(rois[roiNum].t, rois[roiNum].tmin);
					rois[roiNum].w = Math.min(rois[roiNum].w, rois[roiNum].wmax);
					rois[roiNum].h = Math.min(rois[roiNum].h, rois[roiNum].hmax);
					rois[roiNum].t = Math.min(rois[roiNum].t, rois[roiNum].tmax);

					// Repaint
					draw();
				}
			});

			return r;

		}   // end of buildResizeBox()

		// Build translation box
		private Rectangle buildTranslationBox(double x, double y) {

			// Create the rectangle
			Rectangle r = new Rectangle();
			r.setX(x);
			r.setY(y);
			r.setWidth(1.5*cornerSize);
			r.setHeight(1.5*cornerSize);
			r.setFill(Color.BLUE);
			r.setStrokeWidth(1);

			// Make it draggable
			r.addEventHandler(MouseEvent.MOUSE_PRESSED, new EventHandler<MouseEvent>() {
				@Override public void handle(MouseEvent event) {
					setMouse(event.getSceneX(), event.getSceneY()); // Saves (x,y) coordinates of where mouse is clicked in (this.x, this.y)
				}
			});

			r.addEventHandler(MouseEvent.MOUSE_DRAGGED, new EventHandler<MouseEvent>() {
				@Override public void handle(MouseEvent event) {

					// Get the mouse deltas. Account for overall scene scaling (zoom).
					double dx = (event.getSceneX() - getMouseX())/zoom;  // getMouseX() returns this.x
					double dy = (event.getSceneY() - getMouseY())/zoom;

					// Save the current mouse value
					setMouse(event.getSceneX(), event.getSceneY());

					// Update ROI
					if (orientation.equals("Axi")) {
						rois[roiNum].x += dx;
						rois[roiNum].y += dy;
					}
					if (orientation.equals("Sag")) {
						rois[roiNum].y += dx;
						rois[roiNum].z += dy;
					}
					if (orientation.equals("Cor")) {
						rois[roiNum].x += dx;
						rois[roiNum].z += dy;
					}

					/*
					rois[roiNum].x = Math.min(rois[roiNum].x, Meta.nx - rois[roiNum].w/2);
					rois[roiNum].x = Math.max(rois[roiNum].x, rois[roiNum].w/2);
					rois[roiNum].y = Math.min(rois[roiNum].y, Meta.ny - rois[roiNum].h/2);
					rois[roiNum].y = Math.max(rois[roiNum].y, rois[roiNum].h/2);
					rois[roiNum].z = Math.min(rois[roiNum].z, Meta.nz - rois[roiNum].t/2);
					rois[roiNum].z = Math.max(rois[roiNum].z, rois[roiNum].t/2);
					*/

					draw();
				}
			});

			return r;

		}   // end of buildTranslationBox()
	} // End of SelectionBox class

	public static void main(String[] args) {
		launch(args);
	}
}
