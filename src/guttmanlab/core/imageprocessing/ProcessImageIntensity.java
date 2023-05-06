package guttmanlab.core.imageprocessing;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class ProcessImageIntensity {

	public ProcessImageIntensity(File file){
		BufferedImage img = null;
			try {
				img = ImageIO.read(file);
			} catch (IOException e) { e.printStackTrace();}
			double[] val=new double[1];
			
			val=img.getRaster().getPixel(1266, 211, val);
			
			System.err.println(val.length);
		
	}
	
	public static void main(String[] args){
		File file=new File("~/Downloads/Cy5_ZProject_880_11-8-18/1-1_Cy5.czi.tiff");
		new ProcessImageIntensity(file);
	}
	
}
