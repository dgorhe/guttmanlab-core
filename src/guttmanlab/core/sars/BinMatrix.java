package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class BinMatrix {

	public BinMatrix(MatrixWithHeaders mwh, String save, int bin, double minEnrichment) throws IOException{
		
		//bin rows
		//Collection<SingleInterval> regions=getRegions(mwh);
		//Collection<SingleInterval> bins=bin(regions, bin);
		
		MatrixWithHeaders binnedMatrix=binMatrix(mwh, bin, minEnrichment);
		
		
		binnedMatrix.write(save);
	}
	
	
	
	private Collection<SingleInterval> bin(Collection<SingleInterval> regions, int numberOfBins) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: regions){
			
			//split region into number of bins
			int length=region.getLength()/numberOfBins;
			
			for(int i=0; i<numberOfBins; i++){
				int start=i*length;
				int end=start+length;
				SingleInterval bin=new SingleInterval(region.getReferenceName(), start, end);
				rtrn.add(bin);
			}
			
		}
		
		return rtrn;
	}

	/*private Collection<SingleInterval> bin(Collection<SingleInterval> regions, int bin) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: regions){
			Iterator<DerivedAnnotation<? extends Annotation>> iter=region.getGenomicWindows(bin, bin).sortedIterator();
			while(iter.hasNext()){
				DerivedAnnotation a=iter.next();
				SingleInterval i=new SingleInterval(a.getReferenceName(), a.getReferenceStartPosition(), a.getReferenceEndPosition());
				rtrn.add(i);
			}
		}
		
		return rtrn;
	}*/



	private Collection<SingleInterval> getRegions(MatrixWithHeaders mwh) {
		Map<String, Integer> max=new TreeMap<String, Integer>();
		for(String row: mwh.getRowNames()){
			String rna=row.split(":")[0];
			Integer pos=new Integer(row.split(":")[1]);
			if(max.containsKey(rna)){
				int val=max.get(rna);
				pos=Math.max(val, pos);
			}
			max.put(rna, pos);
		}
		
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		for(String rna: max.keySet()){
			int end=max.get(rna);
			SingleInterval region=new SingleInterval(rna, 0, end);
			rtrn.add(region);
		}
		return rtrn;
	}



	private MatrixWithHeaders binMatrix(Collection<SingleInterval> bins, MatrixWithHeaders mwh) {
		List<String> newRows=getRows(bins);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(newRows, mwh.getColumnNames());
		
		for(String column: mwh.getColumnNames()){
			for(SingleInterval bin: bins){
				String row=bin.toUCSC();
				
				List<Double> list=new ArrayList<Double>();
				for(int i=bin.getReferenceStartPosition(); i<bin.getReferenceEndPosition(); i++){
					String oldRow=bin.getReferenceName()+":"+i;
					if(mwh.containsRow(oldRow)){
						double val=mwh.get(oldRow, column);
						list.add(val);
					}
				}
				double p=Statistics.mean(list);
				rtrn.set(row, column, p);
			}
		}
		
		return rtrn;
	}
	
	private MatrixWithHeaders binMatrix(MatrixWithHeaders mwh, int numBins, double minEnrichment) {
		Map<String, Collection<String>> bins=getBins(mwh, numBins);
		
		
		
		bins=filterGenes(mwh, bins, minEnrichment);
		List<String> newRows=getNames(bins);
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(newRows, mwh.getColumnNames());
		
		for(String column: mwh.getColumnNames()){
			for(String newRow: bins.keySet()){
				Collection<String> oldRows=bins.get(newRow);
				List<Double> list=new ArrayList<Double>();
				for(String row: oldRows){
					double val=mwh.get(row, column);
					list.add(val);
				}
				//double p=Statistics.max(list);
				//double p=secondToMax(list);
				double p=Statistics.mean(list);
				//System.err.println(newRow+" "+list.size()+" "+p+" "+oldRows.size());
				rtrn.set(newRow, column, p);
			}
		}
		
		return rtrn;
	}



	private double secondToMax(List<Double> list) {
		Collections.sort(list);
		int index=list.size()-2;
		if(index<0){index=0;}
		return list.get(index);
	}



	private Map<String, Collection<String>> filterGenes(MatrixWithHeaders mwh, Map<String, Collection<String>> bins, double minEnrichment) {
		
		ArrayList<String> list=new ArrayList<String>();
		for(String gene: bins.keySet()){
			double max=0;
			for(String row: bins.get(gene)){
				//String name=gene+":"+row;
				//System.err.println(row);
				double[] vals=mwh.getRow(row);
				max=Math.max(max, Statistics.max(vals));
			}
			if(max>minEnrichment){list.add(gene);}
			//else{System.err.println(gene+" "+max);}
		}
		
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String gene: list){
			rtrn.put(gene, bins.get(gene));
		}
		
		return rtrn;
	}



	private Map<String, Collection<String>> getBins(MatrixWithHeaders mwh, int bins) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		Map<String, TreeSet<Integer>> mapByRNA=getListByRNA(mwh.getRowNames());
		
		for(String rna: mapByRNA.keySet()){
			int numPerBin=(mapByRNA.get(rna).size()/bins);
			if(numPerBin==0){numPerBin=1;}
			
			Iterator<Integer> iter=mapByRNA.get(rna).iterator();
			
			for(int i=0; i<bins; i++){
				Collection<String> binList=new ArrayList<String>();
				while(iter.hasNext()&& binList.size()<numPerBin){
					Integer next=iter.next();
					String name=rna+":"+next;
					binList.add(name);
				}
				String comboName=rna+":"+i;
				if(binList.size()>0){
					rtrn.put(comboName, binList);
				}
			}
			
			
		}
		
		return rtrn;
		
	}



	private Map<String, TreeSet<Integer>> getListByRNA(List<String> rowNames) {
		Map<String, TreeSet<Integer>> rtrn=new TreeMap<String, TreeSet<Integer>>();
		
		for(String row: rowNames){
			String rna=row.split(":")[0];
			Integer pos=new Integer(row.split(":")[1]);
			//System.err.println(row +" "+rna+" "+pos);
			TreeSet<Integer> temp=new TreeSet<Integer>();
			if(rtrn.containsKey(rna)){temp=rtrn.get(rna);}
			temp.add(pos);
			rtrn.put(rna, temp);
		}
		
		return rtrn;
	}



	private List<String> getNames(Map<String, Collection<String>> bins) {
		ArrayList<String> rtrn=new ArrayList<String>();
		rtrn.addAll(bins.keySet());
		return rtrn;
	}



	private List<String> getRows(Collection<SingleInterval> bins) {
		Collection<String> temp=new TreeSet<String>();
		for(SingleInterval bin: bins){temp.add(bin.toUCSC());}
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(temp);
		return rtrn;
	}



	public static void main(String[] args) throws IOException{
		if(args.length>3){
		File file=new File(args[0]);
		String save=args[1];
		int bin=new Integer(args[2]);
		double minEnrichment=new Double(args[3]);
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(file);
		
		new BinMatrix(mwh, save, bin, minEnrichment);
		}
		else{System.err.println(usage);}
	}
	
	static String usage= " args[0]=file \n args[1]=save \n args[2]=number of bins \n args[3]=minEnrichment";
}
