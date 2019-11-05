// Java class for displaying an int[nx][ny] array scaled to [0,255] as a 2D image.
// Keeps track of slice offset so users of this class can calculate 
// the appropriate im[][] array ("slice") to pass to .setPixels().
// 

import javafx.scene.image.ImageView;  
import javafx.scene.image.WritableImage;
import javafx.scene.image.PixelWriter;
import javafx.scene.image.PixelFormat;

import java.nio.ByteBuffer;

public class MRImageView extends ImageView {

	public double sliceOffset = 0;        // offset of this view plane from iso-center (pixels)

	private boolean enableScroll;

	private double minSliceOffset;
	private double maxSliceOffset;

	private WritableImage wr;
	private PixelWriter pw;
	private PixelFormat<ByteBuffer> pixelFormat;
	private byte[] imRgb;

	public MRImageView(int nx, int ny, double minSliceOffset, double maxSliceOffset, boolean enableScroll) { 
		super();

		this.enableScroll = enableScroll;
		this.minSliceOffset = minSliceOffset;
		this.maxSliceOffset = maxSliceOffset;

		wr = new WritableImage(nx,ny);
		this.setImage(wr);

		pw = wr.getPixelWriter();

		pixelFormat = PixelFormat.getByteRgbInstance();

		// change slice offset when scrolling with mouse wheel/pad
		setOnScroll(e -> {
			double delta = e.getDeltaY();

			if (enableScroll) {
				if (delta > 0) {
					sliceOffset += 2.;
				}
				else if (delta < 0) {
					sliceOffset -= 2.;
				}
			}
			
			sliceOffset = Math.max(minSliceOffset, sliceOffset);
			sliceOffset = Math.min(maxSliceOffset, sliceOffset);
		});
	}

	public MRImageView(int nx, int ny, double minSliceOffset, double maxSliceOffset) { 
		this(nx, ny, minSliceOffset, maxSliceOffset, true);
	}

	public MRImageView(int nx, int ny) { 
		this(nx, ny, 0., 0.);
	}

	// update image (pixels)
	public void setPixels(int[][] im) {
		int nx = im.length;
		int ny = im[0].length;
		imRgb = Utils.grayscaleIm2rgb(im, nx, ny);
		pw.setPixels(0, 0, nx, ny, pixelFormat, imRgb, 0, 3*nx);
	}
}
