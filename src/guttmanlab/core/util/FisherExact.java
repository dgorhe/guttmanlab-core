package guttmanlab.core.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class FisherExact {
	
	public FisherExact(File file1, File file2, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Map<SingleInterval, Pair<Integer>> sample1=parse(file1);
		Map<SingleInterval, Pair<Integer>> sample2=parse(file2);
		
		for(SingleInterval region: sample1.keySet()){
			Pair<Integer> pair1=sample1.get(region);
			Pair<Integer> pair2=sample2.get(region);
			if(pair1!=null && pair2!=null){
				//System.err.println(region.toUCSC()+" "+pair1.toString()+" "+pair2.toString());
				double p=computeP(pair1, pair2);
				writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+p+"\t"+pair1.getValue1()+"\t"+pair1.getValue2()+"\t"+pair2.getValue1()+"\t"+pair2.getValue2()+"\n");
			}
			else{System.err.println("Skipped "+region.toUCSC());}
		}
		
		writer.close();
	}
	

	private Map<SingleInterval, Pair<Integer>> parse(File file) throws IOException {
		Map<SingleInterval, Pair<Integer>> rtrn=new TreeMap<SingleInterval, Pair<Integer>>();
		
		Iterator<String> lines=BEDFileIO.loadLines(file.getAbsolutePath()).iterator();
		lines.next();
		
		while(lines.hasNext()){
			String line=lines.next();
			String[] tokens=line.split("\t");
			String chr=tokens[1];
			int start=new Integer(tokens[2]);
			int end=start+1;
			SingleInterval region=new SingleInterval(chr, start, end);
			Pair<Integer> pair=new Pair<Integer>(new Integer(tokens[3]), new Integer(tokens[4]));
			rtrn.put(region, pair);
		}
		
		return rtrn;
	}


	private double computeP(Pair<Integer> pair1, Pair<Integer> pair2) {
		return Statistics.fisherExact(pair1.getValue1(), pair1.getValue2(), pair2.getValue1(), pair2.getValue2());
	}


	public static void main(String[] args) throws IOException{
		File file1=new File(args[0]);
		File file2=new File(args[1]);
		String save=args[2];
		new FisherExact(file1, file2, save);
	}
	
	
	
}
