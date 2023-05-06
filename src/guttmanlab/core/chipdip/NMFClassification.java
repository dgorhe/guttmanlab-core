package guttmanlab.core.chipdip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class NMFClassification {

	public NMFClassification(File input) {
		//run NMF
		//write coefficients
		//categorize regions (bedgraphs)
	}
	
	
	public static void main(String[] args) throws IOException {
		File inputMatrix=new File(args[0]);
		convert(inputMatrix, args[1]);
		
	}


	private static void convert(File input, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			String key=tokens[0]+":"+tokens[1]+"-"+tokens[2];
			writer.write(key);
			for(int i=3; i<tokens.length; i++) {writer.write("\t"+tokens[i]);}
			writer.write("\n");
		}
		reader.close();
		writer.close();
	}
	
}
