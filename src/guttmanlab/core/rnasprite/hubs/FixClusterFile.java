package guttmanlab.core.rnasprite.hubs;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.rnasprite.Cluster;


public class FixClusterFile {

	private static void fix(String input, String output) throws IOException {
		FileWriter writer=new FileWriter(output);
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(input)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			if(tokens.length>2 && tokens.length<=1001) {
				String name=tokens[0];
				Cluster c=new Cluster(name);
				for(int i=1; i<tokens.length; i++) {
					String chr=tokens[i].split(":")[0];
					int start=Integer.parseInt(tokens[i].split(":")[1]);
					int end=start+100;
					SingleInterval r=new SingleInterval(chr, start, end);
					c.addDNARead(r);
				}
				writer.write(c.toString()+"\n");
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		reader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		String input=args[0];
		String output=args[1];
		fix(input, output);
		System.err.println("done");
	}
	
}
