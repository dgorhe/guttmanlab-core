package guttmanlab.core.spidr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import guttmanlab.core.annotation.SingleInterval;

public class MatrixToBedgraphs {

	private static void write(File input, String saveDir) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		
		FileWriter[] writers=null;
		
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(counter==0){
				writers=new FileWriter[tokens.length];
				for(int i=1; i<tokens.length; i++) {
					writers[i]=new FileWriter(saveDir+"/"+tokens[i]+".bedgraph");
				}
			}
			
			else{
				String region=tokens[0];
				SingleInterval interval=new SingleInterval(region);
				for(int i=1; i<tokens.length; i++) {
					double score=Double.parseDouble(tokens[i]);
					writers[i].write(interval.toBedgraph(score)+"\n");
				}
			}
			
			
			
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		reader.close();
		
		for(int i=1; i<writers.length; i++) {writers[i].close();}
		
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length>1) {
		write(new File(args[0]), args[1]);
		} 
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=matrix \n args[1]=save dir";
	
}
