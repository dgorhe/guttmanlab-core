package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class PullBedgraphFromMatrix {

	
	public static void main(String[] args) throws IOException {
		MatrixWithHeaders data=new MatrixWithHeaders(new File(args[0]));
		String saveDir=args[1];
		
		for(String col: data.getColumnNames()) {
			System.err.println(col);
			write(saveDir+"/"+col+".bedgraph", data, col);
		}
		
	}

	private static void write(String string, MatrixWithHeaders data, String col) throws IOException {
		FileWriter writer=new FileWriter(string);
		
		for(String row: data.getRowNames()) {
			SingleInterval region=new SingleInterval(row);
			double val=data.get(row, col);
			writer.write(region.toBedgraph(val)+"\n");
		}
		
		writer.close();
		
	}
	
}
