package guttmanlab.core.chipdip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class CoverageMatrix {
	int min=100;

	public CoverageMatrix(File[] bedgraphs, Map<String, IntervalTree<SingleInterval>> regions, String save, int min) throws IOException {
		this.min=min;
		List<String> columns=getColumns(bedgraphs, regions);
		List<String> rows=getRows(bedgraphs);
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i<bedgraphs.length; i++) {
			String row=bedgraphs[i].getName();
			Map<SingleInterval, Double> map=parse(bedgraphs[i], regions);
			for(SingleInterval r: map.keySet()) {
				String column=r.toUCSC();
				mwh.set(row, column, map.get(r));
			}
		}
		
		mwh=filterRows(mwh);
		mwh.write(save);
	}
	
	
	private MatrixWithHeaders filterRows(MatrixWithHeaders mwh) {
		Collection<String> rows=new TreeSet<String>();
		
		for(String row: mwh.getRowNames()) {
			double[] vals=mwh.getRow(row);
			if(Statistics.sum(vals)>min) {rows.add(row);}
			else {System.err.println("excluded "+row);}
			
		}
		
		return mwh.submatrixByRowNames(rows);
	}


	private List<String> getRows(File[] bedgraphs) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<bedgraphs.length; i++) {
			rtrn.add(bedgraphs[i].getName());
		}
		
		return rtrn;
	}


	private List<String> getColumns(File[] bedgraphs, Map<String, IntervalTree<SingleInterval>> regions) throws NumberFormatException, IOException {
		
		Collection<SingleInterval> set=new TreeSet<SingleInterval>();
		for(int i=0; i<bedgraphs.length; i++) {
			set.addAll(parse(bedgraphs[i], regions).keySet());
		}
		
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval r: set) {
			rtrn.add(r.toUCSC());
		}
		return rtrn;
	}


	private Map<SingleInterval, Double> parse(File file, Map<String, IntervalTree<SingleInterval>> regions) throws NumberFormatException, IOException {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		//System.err.println(file.getName());		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String chr=tokens[0];
			if(!tokens[0].startsWith("chr")) {chr="chr"+tokens[0];}
			SingleInterval r=new SingleInterval(chr, Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			double val=Double.parseDouble(tokens[3]);
			if(overlaps(regions, r)) {
				rtrn.put(r, val);
			}
		}
		reader.close();
		return rtrn;
	}


	private boolean overlaps(Map<String, IntervalTree<SingleInterval>> regions, SingleInterval r) {
		if(regions.containsKey(r.getReferenceName())) {
			IntervalTree<SingleInterval> tree=regions.get(r.getReferenceName());
			if(tree.hasOverlappers(r.getReferenceStartPosition(), r.getReferenceEndPosition())) {return true;}
		}
		return false;
	}


	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			File[] files=new File(args[0]).listFiles();
			Map<String, IntervalTree<SingleInterval>> regions=BEDFileIO.loadSingleIntervalTree(args[1]);
			String save=args[2];
			int min=Integer.parseInt(args[3]);
			new CoverageMatrix(files, regions, save, min);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=bedgraphs (directory) \n args[1]=regions \n args[2]=save";
	
}
