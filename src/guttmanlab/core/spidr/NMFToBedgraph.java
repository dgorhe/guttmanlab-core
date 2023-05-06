package guttmanlab.core.spidr;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class NMFToBedgraph {
	
	private static void write(String save, List<String> lines, int i) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		for(String line: lines) {
			if(counter>0) {
				String[] tokens=line.split("\t");
				SingleInterval interval=new SingleInterval(tokens[0]);
				double score=Double.parseDouble(tokens[i]);
				if(score>0) {
					writer.write(interval.toBedgraph(score)+"\t"+interval.getOrientation()+"\n");
				}
			}
			counter++;
		}
		
		writer.close();
	}
	

	public static void main(String[] args) throws IOException {
		List<String> lines=BEDFileIO.loadLines(args[0]);
		String save=args[1];
		
		String[] header=lines.get(0).split("\t");
		for(int i=1; i<header.length; i++) {
			write(save+"."+header[i]+".bedgraph", lines, i);
		}
		
	}

	
}
