package guttmanlab.core.matrixmanipulation;

import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class ExtractRegionFromMatrix {

	public static MatrixWithHeaders extractSubmatrix(MatrixWithHeaders mwh, SingleInterval region) {
		List<String> rows=new ArrayList<String>();
		
		for(String name: mwh.getRowNames()) {
			SingleInterval r=new SingleInterval(name);
			if(r.overlaps(region)) {rows.add(name);}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, rows) ;
		for(String row: rows) {
			for(String column: rows) {
				rtrn.set(row, column, mwh.get(row, column));
			}
		}
		
		return rtrn;
	}
	
	
	
}
