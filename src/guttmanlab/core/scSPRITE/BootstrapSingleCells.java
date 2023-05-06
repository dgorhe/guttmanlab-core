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
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;

public class BootstrapSingleCells {

	public BootstrapSingleCells(Map<String, MatrixWithHeaders> data1, int perms,  String save, String binToUse) throws IOException, InterruptedException{
		//Make matrix
		MatrixWithHeaders observed=contactMatrix(data1);
		
		System.err.println("Constructed observed");
		
		File observedOutput=writeAndNorm(observed, save+".observed.matrix");
		//observed.write();
		observed=new MatrixWithHeaders(observedOutput);
		
		System.err.println("ICE normalized");
		
		List<String> rows=new ArrayList<String>();
		rows.add("observed");
		for(int i=0; i<perms; i++){rows.add("perm"+i);}
		MatrixWithHeaders merged=new MatrixWithHeaders(rows, observed.getColumnNames());
		merged.setRow("observed", observed.getRow(binToUse));
		
		
		
		//MatrixWithHeaders[] iced=new MatrixWithHeaders[perms];
		//File[] icedFiles=new File[perms];
		for(int i=0; i<perms; i++){
			System.err.println("perm"+i);
			MatrixWithHeaders random=permute(data1, observed);
			String input=save+".perm"+i+".matrix";
			File output=writeAndNorm(random, input);
			//icedFiles[i]=output;
			MatrixWithHeaders iced=new MatrixWithHeaders(output);
			System.err.println("ICE norm");
			merged.setRow("perm"+i, iced.getRow(binToUse));
		}
		
		//MatrixWithHeaders merged=mergedMatrix(observed, iced, binToUse);
		merged.write(save+".merged.matrix");
		
	}

	/*private MatrixWithHeaders mergedMatrix(MatrixWithHeaders observed, MatrixWithHeaders[] iced, String binToUse) {
		List<String> rows=new ArrayList<String>();
		rows.add("observed");
		for(int i=0; i<iced.length; i++){rows.add("perm"+i);}
		List<String> columns=observed.getColumnNames();
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		
		
		for(int i=0; i<iced.length; i++){
			rtrn.setRow("perm"+i, iced[i].getRow(binToUse));
		}
		
		return rtrn;
	}*/

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
			if(i%10==0){System.err.println(i+" "+listFiles.length);}
			rtrn.put(listFiles[i].getName(), new MatrixWithHeaders(listFiles[i]));
		}
		System.err.println("Done loading matrices");
		
		return rtrn;
	}

	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
			Map<String, MatrixWithHeaders> data1=parse(new File(args[0]).listFiles());
			int numPerms=new Integer(args[1]);
			String save=args[2];
			String regionToUse=args[3];
			new BootstrapSingleCells(data1, numPerms, save, regionToUse);
		}
		else{System.err.println(usage);}
	}
	
	
	

	static String usage=" args[0]=directory of cluster files \n args[1]=num perms \n args[2]=save \n args[3]=bin to use";
	
}
