package guttmanlab.core.chip;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class MakeTable {

	public MakeTable(File[] bedFiles, int binSize, String save) throws IOException {
		
		List<String> columns=getColumns(bedFiles);
		List<String> rows=getRows(bedFiles, binSize);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		
		for(int i=0; i<bedFiles.length; i++) {
			Map<SingleInterval, Double> vals=parse(bedFiles[i], binSize);
			add(mwh, vals, bedFiles[i].getName());
		}
		
		//mwh=filter(mwh, 2);
		mwh=filterColumns(mwh);
		
		mwh.write(save);
	}

	private MatrixWithHeaders filterColumns(MatrixWithHeaders mwh) {
		List<String> list=new ArrayList<String>();
		for(String column: mwh.getColumnNames()) {
			double[] vals=mwh.getColumn(column);
			if(!filter(vals)) {list.add(column);}
		}
		return mwh.submatrixByColumnNames(list);
	}

	private boolean filter(double[] vals) {
		int count=number(vals);
		return count<100;
	}

	private MatrixWithHeaders filter(MatrixWithHeaders mwh, int number) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()) {
			double[] vals=mwh.getRow(row);
			if(number(vals)>=number) {list.add(row);}
		}
		return mwh.submatrixByRowNames(list);
	}

	private int number(double[] vals) {
		int count=0;
		for(int i=0; i<vals.length; i++) {
			if(vals[i]>2.0) {count++;}
		}
		return count;
	}

	private List<String> getRows(File[] bedFiles, int binSize) throws IOException {
		Collection<SingleInterval> bins=new TreeSet<SingleInterval>();
		
		for(int i=0; i<bedFiles.length; i++) {
			Map<SingleInterval, Double> vals=parse(bedFiles[i], binSize);
			bins.addAll(vals.keySet());
		}
		
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval r: bins) {rtrn.add(r.toUCSC());}
		
		return rtrn;
	}

	private List<String> getColumns(File[] bedFiles) {
		List<String> rtrn=new ArrayList<String>();
		for(int i=0; i<bedFiles.length; i++) {rtrn.add(bedFiles[i].getName());}
		return rtrn;
	}

	private Map<SingleInterval, Double> parse(File file, int binSize) throws IOException {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		List<String> lines=BEDFileIO.loadLines(file.getAbsolutePath());
		for(String line: lines) {
			String[] tokens=line.split("\t");
			SingleInterval r=new SingleInterval(tokens[0], Integer.parseInt(tokens[1]), Integer.parseInt(tokens[2]));
			Collection<SingleInterval> bins=r.getWindowsCollection(binSize, binSize);
			double val=Double.parseDouble(tokens[4]);
			for(SingleInterval bin: bins) {
				if(rtrn.containsKey(bin)) {val=Math.max(rtrn.get(bin), val);}
				rtrn.put(bin, val);
			}
		}
		
		return rtrn;
	}

	private void add(MatrixWithHeaders mwh, Map<SingleInterval, Double> vals, String name) {
		for(SingleInterval r: vals.keySet()) {
			String row=r.toUCSC();
			double val=vals.get(r);
			mwh.set(row, name, val);
		}
		
	}
	
	
	public static void main(String[] args) throws IOException {
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		int binSize=100;
		new MakeTable(files, binSize, save);
		
	}
	
}
