package guttmanlab.core.proteinSPRITE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.rnasprite.Cluster;

public class ProteinHeatmap {

	public ProteinHeatmap(BarcodingDataStreaming data, String save) throws IOException, InterruptedException {
		List<String> rows=getProteins(data);
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, rows);
		
		int counter=0;
		int pairs=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			Collection<String> proteins=c.getProteinSet();
			if(proteins.size()>=2) {
				addToMatrix(mwh, proteins);
				pairs++;
			}
			counter++;
			if(counter%1000000==0) {System.err.println(counter+" "+pairs);}
		}
		
		//Map<String, Double> counts=getRowCounts(mwh);
		//mwh=norm(mwh, counts);
		writeAndNorm(mwh, save);
		//mwh.write(save);
		data.close();
	}
	
	
	private static File writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return rtrn;
	}
	

	private MatrixWithHeaders norm(MatrixWithHeaders mwh, Map<String, Double> counts) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		for(String row: mwh.getRowNames()) {
			for(String column: mwh.getColumnNames()) {
				double score=mwh.get(row, column);
				double expected=((counts.get(row)/(double)counts.size())+(counts.get(column)/(double)counts.size()));
				rtrn.set(row, column, score/expected);
			}
		}
		return rtrn;
	}


	private Map<String, Double> getRowCounts(MatrixWithHeaders mwh) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String row: mwh.getRowNames()) {
			double sum=Statistics.sum(mwh.getRow(row));
			rtrn.put(row, sum);
		}
		
		return rtrn;
	}


	private void addToMatrix(MatrixWithHeaders mwh, Collection<String> proteins) {
		for(String p1: proteins) {
			for(String p2: proteins) {
				//if(!p1.equals(p2)) {
					mwh.incrementCount(p1, p2);
				//}
			}
		}
		
	}

	private List<String> getProteins(BarcodingDataStreaming data) {
		Collection<String> proteins=new TreeSet<String>();
		int counter=0;
		while(data.hasNext()) {
			Cluster c=data.next();
			proteins.addAll(c.getProteinSet());
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		data.close();
		
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(proteins);
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new ProteinHeatmap(data, save);
	}
	
}
