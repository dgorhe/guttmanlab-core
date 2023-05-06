package guttmanlab.core.spidr;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class BinMatrix {
	
	private static MatrixWithHeaders bin(MatrixWithHeaders matrix, int binResolution) {
		
		List<String> newRows=binRows(matrix.getRowNames(), binResolution);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(newRows, matrix.getColumnNames());
		
		for(String column: matrix.getColumnNames()) {
			for(String row: matrix.getRowNames()) {
				String binned=bin(row, binResolution);
				double score=matrix.get(row, column);
				double newScore=Math.max(score, rtrn.get(binned, column));
				rtrn.set(binned, column, newScore);
			}
		}
		return rtrn;
	}
	
	

	private static List<String> binRows(List<String> rowNames, int binResolution) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String row: rowNames) {
			String newRow=bin(row, binResolution);
			if(!rtrn.contains(newRow)) {
				rtrn.add(newRow);
			}
			
		}
		
		return rtrn;
	}



	private static String bin(String row, int binResolution) {
		SingleInterval interval=new SingleInterval(row);
		SingleInterval binned=interval.bin(binResolution);
		binned.setOrientation(interval.getOrientation());
		return binned.toUCSCStrand();
	}



	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			MatrixWithHeaders matrix=new MatrixWithHeaders(new File(args[0]), ",");
			int binResolution=Integer.parseInt(args[1]);
			String save=args[2];
			MatrixWithHeaders binned=bin(matrix, binResolution);
			binned.write(save);
		}
		else {System.err.println(usage);}
	}

	

	static String usage=" args[0]=matrix \n args[1]=bin size \n args[2]=save";
	
}
