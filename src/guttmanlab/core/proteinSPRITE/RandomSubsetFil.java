package guttmanlab.core.proteinSPRITE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class RandomSubsetFil {

	public RandomSubsetFil(File file, double fraction, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		
		int count=0;
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			double rand=Math.random();
			if(rand<fraction) {writer.write(nextLine+"\n");}
			
			count++;
			if(count%100000==0) {System.err.println(count);}
		}
		reader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		File f=new File(args[0]);
		double fract=Double.parseDouble(args[1]);
		String save=args[2];
		new RandomSubsetFil(f, fract, save);
	}
	
}
