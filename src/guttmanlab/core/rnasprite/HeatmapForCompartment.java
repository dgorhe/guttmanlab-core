package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.barcoding.analysis.Cluster;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.simulation.CoordinateSpace;

public class HeatmapForCompartment {

	int binSize=20;
	
	public HeatmapForCompartment(MatrixWithHeaders mwh, SingleInterval region, int binResolution, String save, int binSize) throws IOException {
		this.binSize=binSize;
		
		Collection<SingleInterval> bins=getBins(region, binResolution);
		List<String> rows=getString(bins);
		double normFactor=getIntrachromosomalAverage(mwh);
		double normFactor2=getWindowAverage(mwh, binSize);
		System.err.println(normFactor+" "+normFactor2);
		
		
		MatrixWithHeaders subMatrix=new MatrixWithHeaders(rows, rows);
		
		for(String r1: rows){
			for(String r2: rows){
				if(!r1.equals(r2)) {
					if(mwh.containsColumn(r2) && mwh.containsRow(r1)) {
						subMatrix.set(r1, r2, mwh.get(r1, r2));
					}
				}
			}
		}
		
		subMatrix=normalize(subMatrix, normFactor2);
		
		subMatrix.write(save);
		
	}
	
	public static MatrixWithHeaders norm(MatrixWithHeaders mwh, int binResolution) throws IOException {
		double normFactor2=getWindowAverage(mwh, 10);
		
		
		MatrixWithHeaders subMatrix=new MatrixWithHeaders(mwh.getRowNames(), mwh.getColumnNames());
		
		for(String row: mwh.getRowNames()){
			for(String column: mwh.getColumnNames()){
				subMatrix.set(row, column, mwh.get(row, column)/normFactor2);
					
			}
		}
		
		
		return subMatrix;
	}
	
	public static MatrixWithHeaders norm(MatrixWithHeaders mwh, int binResolution, SingleInterval region) throws IOException {
		double normFactor2=getWindowAverage(mwh, 10);
		
		
		Collection<SingleInterval> bins=getRegions(region, binResolution);
		List<String> rows=getString(bins);
		
		
		
		MatrixWithHeaders subMatrix=new MatrixWithHeaders(rows, rows);
		
		for(String r1: rows){
			for(String r2: rows){
				if(!r1.equals(r2)) {
					if(mwh.containsColumn(r2) && mwh.containsRow(r1)) {
						subMatrix.set(r1, r2, mwh.get(r1, r2)/normFactor2);
					}
				}
			}
		}
		
		return subMatrix;
	}
	
	private MatrixWithHeaders normalize(MatrixWithHeaders subMatrix, double normFactor) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(subMatrix.getRowNames(), subMatrix.getColumnNames());
		
		for(String row: subMatrix.getRowNames()) {
			for(String column: subMatrix.getColumnNames()) {
				rtrn.set(row, column, subMatrix.get(row, column)/normFactor);
			}
		}
		
		return rtrn;
	}

	
	private static double getWindowAverage(MatrixWithHeaders matrix, int binSize) {
		//get window average
		double sum=0;
		double count=0;
		
		
		Iterator<Collection<String>> iter=rowIterator(binSize, matrix);
		
		while(iter.hasNext()) {
			Collection<String> rowsAndCols=iter.next();
			for(String row: rowsAndCols) {
				for(String column: rowsAndCols) {
					if(!row.equals(column)) {
						sum+=matrix.get(row, column);
						count++;
					}
				}
			}
		}
		
		double average=sum/count;
		return average;
	}
	
	private static Iterator<Collection<String>> rowIterator(int binSize, MatrixWithHeaders matrix) {
		Map<String, List<String>> listByChromosome=splitByChromosome(matrix.getRowNames());
		
		List<Collection<String>> rtrn=new ArrayList<Collection<String>>();
		
		for(String chr: listByChromosome.keySet()) {
			System.err.println(chr+" "+listByChromosome.get(chr).size());
			List<String> regions=listByChromosome.get(chr);
			for(int i=0; i<regions.size(); i+=binSize) {
				Collection<String> set=makeCollection(regions, i, i+binSize);
				rtrn.add(set);
			}
		}
		return rtrn.iterator();
	}

	private static Map<String, List<String>> splitByChromosome(List<String> rowNames) {
		Map<String, List<String>> rtrn=new TreeMap<String, List<String>>();
		
		for(String row: rowNames) {
			String chr=row.split(":")[0];
			if(!rtrn.containsKey(chr)) {rtrn.put(chr, new ArrayList<String>());}
			List<String> list=rtrn.get(chr);
			list.add(row);
		}
		
		return rtrn;
	}

	private static Collection<String> makeCollection(List<String> regions, int start, int end) {
		Collection<String> rtrn=new TreeSet<String>();
		
		end=Math.min(end, regions.size());
		
		for(int i=start; i<end; i++) {rtrn.add(regions.get(i));}
		
		return rtrn;
	}

	private double getIntrachromosomalAverage(MatrixWithHeaders matrix) {
		//get average intra-chromosomal
		double sum=0;
		double count=0;
		
		for(String row: matrix.getRowNames()) {
			for(String column: matrix.getColumnNames()) {
				if(!row.equals(column)) {
					if(row.split(":")[0].equals(column.split(":")[0])) {
						sum+=matrix.get(row, column);
						count++;
					}
				}
			}
		}
		return sum/count;
	}

	private static List<String> getString(Collection<SingleInterval> bins) {
		List<String> rtrn=new ArrayList<String>();
		
		for(SingleInterval bin: bins) {rtrn.add(bin.toUCSC());}
		
		return rtrn;
	}

	private Collection<Cluster> generatePairs(Collection<SingleInterval> regions) {
		Collection<Cluster> rtrn=new TreeSet<Cluster>();
		for(SingleInterval region1: regions){
			for(SingleInterval region2: regions){
				if(!region1.getReferenceName().equals(region2.getReferenceName())){
					Cluster c=new Cluster("pairs");
					c.addRead(region1);
					c.addRead(region2);
					rtrn.add(c);
				}
			}
		}
		return rtrn;
	}
	
	
	private Collection<SingleInterval> getBins(Collection<SingleInterval> regionsCombined, int binResolution) {
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		for(SingleInterval region: regionsCombined){
			rtrn.addAll(getRegions(region, binResolution));
		}
		
		return rtrn;
	}
	
	private Collection<SingleInterval> getBins(SingleInterval region, int binResolution) {
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		rtrn.addAll(getRegions(region, binResolution));
		
		
		return rtrn;
	}
	
	private static Collection<SingleInterval> getRegions(SingleInterval region, int binResolution) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(int i=region.getReferenceStartPosition(); i<region.getReferenceEndPosition(); i+=binResolution){
			SingleInterval temp=new SingleInterval(region.getReferenceName(), i, i+binResolution);
			rtrn.add(temp);
			//System.err.println(region1.toUCSC()+" "+region.toUCSC());
		}
		return rtrn;
	}
	
	

	public static void main(String[] args) throws IOException {
		if(args.length>4) {
			MatrixWithHeaders mwh=new MatrixWithHeaders(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			int binResolution=Integer.parseInt(args[2]);
			String save=args[3];
			int binSize=Integer.parseInt(args[4]);
			new HeatmapForCompartment(mwh, region, binResolution, save, binSize);
		}
		else {System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=matrix \n args[1]=regions \n args[2]=bin resolution \n args[3]=save \n args[4]=bin size (for normalization)";
	
}
