package guttmanlab.core.proteinSPRITE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class CountsByCluster {

	
	private static void filter(File in, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(in)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			int size=tokens.length-1;
			writer.write(tokens[0]+"\t"+size+"\n");
		}
		reader.close();
		writer.close();
		
		
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>1) {
			File in=new File(args[0]);
			
			String save=args[1];
			
			filter(in, save);
		}
		else {System.err.println(usage);}
		
	}
	
	

	static String usage=" args[0]=clusters \n args[1]=save";
}
