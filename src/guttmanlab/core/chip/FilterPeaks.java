package guttmanlab.core.chip;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class FilterPeaks {

	public static void main(String[] args) throws NumberFormatException, IOException {
		double threshold=0;
		if(args.length>3) {threshold=Double.parseDouble(args[3]);}
		filterByLength(args[0], args[1], Integer.parseInt(args[2]), threshold);
	}

	private static void filterByLength(String input, String output, int length, double threshold) throws IOException {
		FileWriter writer=new FileWriter(output);
		List<String> lines=BEDFileIO.loadLines(input);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			SingleInterval region=new SingleInterval(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			double score=Double.parseDouble(tokens[3]);
			if(region.size()<=length && score>threshold) {writer.write(line+"\n");}
		}
		writer.close();
	}
	
	
}
