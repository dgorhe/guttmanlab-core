package guttmanlab.core.sars;

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
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class MakeMatrix {

	public MakeMatrix(File[] files, String save) throws IOException{
		Map<String, Double>[] maps=new Map[files.length];
		for(int i=0; i<files.length; i++){
			System.err.println(files[i]);
			maps[i]=parse(files[i]);
		}
		
		MatrixWithHeaders mwh=makeMatrix(maps, files);
		
		//mwh=filter(mwh, val);
		
		mwh.write(save);
	}
	
	private MatrixWithHeaders filter(MatrixWithHeaders mwh, double d) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()){
			double[] vals=mwh.getRow(row);
			if(Statistics.max(vals)>d){list.add(row);}
		}
		return mwh.submatrixByRowNames(list);
	}

	private Map<String, Double> parse(File file) throws NumberFormatException, IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			String key=tokens[0]+":"+tokens[1];
			Double val=new Double(tokens[3]);
			if(Double.isInfinite(val)){val=0.0;}
			rtrn.put(key, val);
		}
		reader.close();
		return rtrn;
	}

	private MatrixWithHeaders makeMatrix(Map<String, Double>[] maps, File[] files) {
		Collection<String> set=new TreeSet<String>();
		List<String> columns=new ArrayList<String>();
		for(int i=0; i< files.length; i++){
			columns.add(files[i].getName());
			set.addAll(maps[i].keySet());
		}
		List<String> rows=new ArrayList<String>();
		rows.addAll(set);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i<maps.length; i++){
			String column=files[i].getName();
			for(String pos: maps[i].keySet()){
				double val=maps[i].get(pos);
				mwh.set(pos, column, val);
			}
		}
		
		return mwh;
	}

	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		new MakeMatrix(files, save);
	}
	
}
