package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class MergeMatrices {

	public static void main (String[] args) throws IOException{
		MatrixWithHeaders m1=new MatrixWithHeaders(new File(args[0]));
		MatrixWithHeaders m2=new MatrixWithHeaders(new File(args[1]));
		
		MatrixWithHeaders merged=merge(m1, m2);
		merged=filter(merged);
		merged.write(args[2]);
	}

	private static MatrixWithHeaders filter(MatrixWithHeaders merged) {
		Collection<String> list=new ArrayList<String>();
		for(String column: merged.getColumnNames()){
			double[] vals=merged.getColumn(column);
			if(Statistics.sum(vals)>0){list.add(column);}
		}
		return merged.submatrixByColumnNames(list);
	}

	private static MatrixWithHeaders merge(MatrixWithHeaders m1, MatrixWithHeaders m2) {
		Collection<String> columns=new TreeSet<String>();
		columns.addAll(m1.getColumnNames());
		columns.addAll(m2.getColumnNames());
		
		List<String> columnNames=new ArrayList<String>();
		columnNames.addAll(columns);
		
		List<String> rows=new ArrayList<String>();
		rows.addAll(m1.getRowNames());
		rows.addAll(m2.getRowNames());
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columnNames);
		
		for(String row: m1.getRowNames()){
			for(String column: m1.getColumnNames()){
				double val=m1.get(row, column);
				rtrn.set(row, column, val);
			}
		}
		
		for(String row: m2.getRowNames()){
			for(String column: m2.getColumnNames()){
				double val=m2.get(row, column);
				rtrn.set(row, column, val);
			}
		}
		
		return rtrn;
	}
	
}
