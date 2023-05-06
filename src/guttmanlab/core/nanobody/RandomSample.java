package guttmanlab.core.nanobody;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;

public class RandomSample {

	public static void main(String[] args) throws IOException {
		File input=new File(args[0]);
		double fraction=Double.parseDouble(args[1]);
		String save=args[2];
		downsample(input, fraction, save);
		
	}

	private static void downsample(File input, double fraction, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			double rand=Math.random();
			if(rand<fraction) {
				writer.write(nextLine+"\n");
			}
			
			
			if(counter%1000000==0) {System.err.println(counter);}
			counter++;
		}
		reader.close();
		writer.close();
	}
	
}
