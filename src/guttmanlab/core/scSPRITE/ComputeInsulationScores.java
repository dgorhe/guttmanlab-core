package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class ComputeInsulationScores {

	int windowSize=5;
	
	public ComputeInsulationScores(MatrixWithHeaders matrix, String save) throws IOException {
	
		FileWriter writer=new FileWriter(save);
		
		//“–ss 80000 –im iqrMean –is 480000 –ids 320000” with contact maps binned at 40kb resolution
		
		for(String name: matrix.getRowNames()) {
			SingleInterval region=new SingleInterval(name);
			List<String> A=getLeft(matrix, name);
			List<String> B=getRight(matrix, name);
			double insulationScore=getScore(matrix, A, B);
			writer.write(region.toBedgraph(insulationScore)+"\n");
		}
		
		writer.close();
	}

	private List<String> getRight(MatrixWithHeaders matrix, String name) {
		int index=matrix.getRowNames().indexOf(name);
			
		List<String> rtrn=new ArrayList<String>();
		for(int i=1; i<windowSize; i++) {
			if(index+i<matrix.getRowNames().size()) {
				String a=matrix.getRowName(index+i);
				rtrn.add(a);
			}
		}
			
		return rtrn;
	}
	

	private List<String> getLeft(MatrixWithHeaders matrix, String name) {
		int index=matrix.getRowNames().indexOf(name);
		
		List<String> rtrn=new ArrayList<String>();
		for(int i=1; i<windowSize; i++) {
			if(index-i>=0) {
				String a=matrix.getRowName(index-i);
				rtrn.add(a);
			}
		}
		
		return rtrn;
	}

	private double getScore(MatrixWithHeaders matrix, List<String> a, List<String> b) {
		List<Double> vals=new ArrayList<Double>();
		for(String a1: a) {
			for(String b1: b) {
				double val=matrix.get(a1, b1);
				vals.add(val);
			}
		}
		return Statistics.mean(vals);
	}
	
	private static MatrixWithHeaders parseHiC(String file, String chr) throws IOException {
		
		List<String> lines=BEDFileIO.loadLines(file);
		List<String> rows=makeRows(lines, chr);
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, rows);
				
		for(int rowIndex=0; rowIndex< lines.size(); rowIndex++){
			String[] numbers=lines.get(rowIndex).trim().split("\t");
			//System.err.println(numbers.length);
			String row=rowIndex+"|mm9|"+chr+":"+(rowIndex*40000)+"-"+((rowIndex+1)*40000);
			for(int columnIndex=0; columnIndex<numbers.length; columnIndex++){
					String column=columnIndex+"|mm9|"+chr+":"+(columnIndex*40000)+"-"+((columnIndex+1)*40000);
					//System.err.println(row +" "+column+" "+columnIndex+" "+numbers[columnIndex]);
					mwh.set(row, column, Double.parseDouble(numbers[columnIndex]));
				}
			}
				
			
			
			return mwh;
		}

	
	private static File writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return rtrn;
	}
		
		
		private static List<String> makeRows(List<String> lines, String chr) {
			List<String> rtrn=new ArrayList<String>();
			for(int i=0; i<lines.size(); i++){
				String row=i+"|mm9|"+chr+":"+(i*40000)+"-"+((i+1)*40000);
				rtrn.add(row);
			}
			return rtrn;
		}
		
		
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>2) {
		MatrixWithHeaders mwh=parseHiC(args[0], args[2]);
		String save=args[1];
		writeAndNorm(mwh, save);
		}else {System.err.println(usage);}
	}

	static String usage=" args[0]=raw matrix \n args[1]=save \n args[2]=chr";
	
	
}
