package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.datastructures.MatrixWithHeaders;

public class SplitSymetricMatrix {

	public SplitSymetricMatrix(MatrixWithHeaders matrix1, MatrixWithHeaders matrix2, String save) throws IOException {
		//double max1=max(matrix1);
		//double max2=max(matrix2);
		
		double max1=1;
		double max2=1;
		
		MatrixWithHeaders merged=new MatrixWithHeaders(matrix1.getRowNames(), matrix1.getColumnNames());
		
		
		for(String row: matrix1.getRowNames()) {
			boolean done=false;
			for(String column: matrix1.getColumnNames()) {
				if(row.equals(column)) {done=true;}
				if(!done) {merged.set(row, column, matrix1.get(row, column)/max1);}
			}
		}
		
		
		for(String row: matrix2.getRowNames()) {
			boolean start=false;
			for(String column: matrix2.getColumnNames()) {
				if(row.equals(column)) {start=true;}
				if(start) {
					if(merged.containsRow(row) && merged.containsColumn(column)) {
					merged.set(row, column, matrix2.get(row, column)/max2);}
					}
			}
		}
		
		merged.write(save);
		
	}
	
	
	private double max(MatrixWithHeaders matrix1) {
		double max=0;
		for(String row: matrix1.getRowNames()) {
			for(String column: matrix1.getColumnNames()) {
				max=Math.max(matrix1.get(row, column), max);
			}
		}
		return max;
	}


	public static void main(String[] args) throws IOException {
		MatrixWithHeaders m1=new MatrixWithHeaders(new File(args[0]));
		MatrixWithHeaders m2=new MatrixWithHeaders(new File(args[1]));
		String save=args[2];
		new SplitSymetricMatrix(m1, m2, save);
	}
	
}
