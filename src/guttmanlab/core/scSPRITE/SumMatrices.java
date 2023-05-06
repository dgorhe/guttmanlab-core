package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;

public class SumMatrices {

	public SumMatrices(Map<String, MatrixWithHeaders> data1, String save) throws IOException, InterruptedException{
		//Make matrix
		MatrixWithHeaders observed=contactMatrix(data1);
		
		System.err.println("Constructed observed");
		
		writeAndNorm(observed, save+".observed.matrix");
		
	}

	

	private File writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return rtrn;
	}

	private MatrixWithHeaders permute(Map<String, MatrixWithHeaders> data1, MatrixWithHeaders observed) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		
		ArrayList<String> names=new ArrayList<String>();
		names.addAll(data1.keySet());
		
		for(int i=0; i<names.size(); i++){
			int index=new Double(Math.random()*names.size()).intValue();
			String name=names.get(index);
			MatrixWithHeaders data=data1.get(name);
			merge(data, rtrn);
		}
		
		
		return rtrn;
		
		
		
	}
	
	private ArrayList<BarcodingDataStreaming> permute(ArrayList<BarcodingDataStreaming> data1) {
		ArrayList<BarcodingDataStreaming> perm=new ArrayList<BarcodingDataStreaming>();
		
		for(int i=0; i<data1.size(); i++){
			int index=new Double(Math.random()*data1.size()).intValue();
			perm.add(data1.get(index));
		}
		
		return perm;
	}
	
	private MatrixWithHeaders contactMatrix(Map<String, MatrixWithHeaders> data1) {
		List<String> rows=getRows(data1);
		List<String> columns=getColumns(data1);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(String name: data1.keySet()){
			MatrixWithHeaders data=data1.get(name);
			merge(data, rtrn);
		}
		
		return rtrn;
	}
	
	
	private List<String> getColumns(Map<String, MatrixWithHeaders> data1) {
		Set<String> rtrn=new TreeSet<String>();
		for(String name: data1.keySet()){
			rtrn.addAll(data1.get(name).getColumnNames());
		}
		
		List<String> list=new ArrayList<String>();
		list.addAll(rtrn);
		return list;
	}

	private List<String> getRows(Map<String, MatrixWithHeaders> data1) {
		Set<String> rtrn=new TreeSet<String>();
		for(String name: data1.keySet()){
			rtrn.addAll(data1.get(name).getRowNames());
		}
		
		List<String> list=new ArrayList<String>();
		list.addAll(rtrn);
		return list;
	}

	private void merge(MatrixWithHeaders data, MatrixWithHeaders rtrn) {
		for(String row: data.getRowNames()){
			for(String column: data.getColumnNames()){
				double score1=rtrn.get(row, column);
				double score2=data.get(row, column);
				double sum=score1+score2;
				rtrn.set(row, column, sum);
			}
		}
		
	}

	private MatrixWithHeaders contactMatrix(Collection<BarcodingDataStreaming> data1, int binResolution, boolean weight, SingleInterval region) {
		MatrixWithHeaders rtrn=getDNADNAContactMatrix(data1, binResolution, weight, region);
		return rtrn;
	}
	
	
	public MatrixWithHeaders getDNADNAContactMatrix(Collection<BarcodingDataStreaming> data1, int binResolution, boolean weight, SingleInterval region){
		List<String> regions=getGenomePositions(binResolution, region);
		System.err.println(regions);
		
		MatrixWithHeaders counts=new MatrixWithHeaders(regions, regions);
		
		for(BarcodingDataStreaming data: data1){
			System.err.println(data.getBarcodeFile().getName());
			scoreDNADNA(data, counts, binResolution, weight);
		}
		return counts;
	}
	
	
	private void scoreDNADNA(BarcodingDataStreaming data, MatrixWithHeaders counts, int binResolution, boolean weight) {
		Collection<Cluster> clusters=data.getClusters();
		for(Cluster c: clusters){
			Cluster binned=c.bin(binResolution);
			for(SingleInterval region1: binned.getAllDNAIntervals()){
				for(SingleInterval region2: binned.getAllDNAIntervals()){
					String row=region1.toUCSC();
					String column=region2.toUCSC();
					if(counts.containsColumn(column) && counts.containsRow(row)){
						//System.err.println(row+" "+column);
						double count=counts.get(row, column);
						double score=2.0/c.getAllDNAIntervals().size();
						//count++;
						
						if(weight){
							count+=score;
						}
						else{
							count++;
						}
						//System.err.println(count);
						counts.set(row, column, count);
					}
				}
			}
		
		}
		
		
		
	}
	
	
	
	private List<String> getGenomePositions(int binResolution, SingleInterval region) {
		ArrayList<String> list=new ArrayList<String>();
		
		SingleInterval start=region.bin(binResolution);
		SingleInterval end=new SingleInterval(region.getReferenceName(), region.getReferenceEndPosition()-1, region.getReferenceEndPosition()).bin(binResolution);
		
		for(int i=start.getReferenceStartPosition(); i<end.getReferenceEndPosition(); i+=binResolution){
			SingleInterval bin=new SingleInterval(start.getReferenceName(), i, i+binResolution);
			list.add(bin.toUCSC());
		}
		
		return list;
	}
	
	
	private static ArrayList<BarcodingDataStreaming> parse(File[] listFiles, SingleInterval region) throws IOException {
		ArrayList<BarcodingDataStreaming> rtrn=new ArrayList<BarcodingDataStreaming>();
		for(int i=0; i<listFiles.length; i++){
			rtrn.add(new BarcodingDataStreaming(listFiles[i], region));
		}
		return rtrn;
	}

	
	private static Map<String, MatrixWithHeaders> parse(File[] listFiles) throws IOException {
		Map<String, MatrixWithHeaders> rtrn=new TreeMap<String, MatrixWithHeaders>();
		
		for(int i=0; i<listFiles.length; i++){
			if(listFiles[i].getName().endsWith(".weighted")){
			if(i%10==0){System.err.println(i+" "+listFiles.length);}
			
			List<String> lines=BEDFileIO.loadLines(listFiles[i].getAbsolutePath());
			List<String> rows=makeRows(lines);
			MatrixWithHeaders mwh=new MatrixWithHeaders(rows, rows);
			
			for(int rowIndex=0; rowIndex< lines.size(); rowIndex++){
				String[] numbers=lines.get(rowIndex).split("\t");
				String row="chr6:"+(rowIndex*40000)+"-"+((rowIndex+1)*40000);
				for(int columnIndex=0; columnIndex<numbers.length; columnIndex++){
					String column="chr6:"+(columnIndex*40000)+"-"+((columnIndex+1)*40000);
					mwh.set(row, column, new Double(numbers[columnIndex]));
				}
			}
			
			rtrn.put(listFiles[i].getName(), mwh);
			}
		}
		System.err.println("Done loading matrices");
		
		return rtrn;
	}

	
	
	private static List<String> makeRows(List<String> lines) {
		List<String> rtrn=new ArrayList<String>();
		for(int i=0; i<lines.size(); i++){
			String row="chr6:"+(i*40000)+"-"+((i+1)*40000);
			rtrn.add(row);
		}
		return rtrn;
	}



	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>1){
			Map<String, MatrixWithHeaders> data1=parse(new File(args[0]).listFiles());
			String save=args[1];
			new SumMatrices(data1, save);
		}
		else{System.err.println(usage);}
	}
	
	
	

	static String usage=" args[0]=directory of cluster files \n args[1]=save";
	
}
