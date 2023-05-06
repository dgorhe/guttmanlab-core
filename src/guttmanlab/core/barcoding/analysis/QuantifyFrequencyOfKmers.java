package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class QuantifyFrequencyOfKmers {

	public QuantifyFrequencyOfKmers(BarcodingData data, Collection<Cluster> clusters) throws NumberFormatException, IOException{
		for(Cluster cluster: clusters){
			int count=quantify(cluster, data);
			System.out.println(cluster.toStringNoName()+"\t"+count);
		}
	}

	public static int quantify(Cluster cluster, BarcodingData data) throws NumberFormatException, IOException {
		Collection<String> merged=null;
		for(SingleInterval region: cluster.getAllIntervals()){
			Collection<String> barcodes=data.getBarcodesOverlappingRegion(region);
			merged=merge(merged, barcodes);
		}
		return merged.size();
	}
	
	public static Collection<String> getClusters(Cluster cluster, BarcodingData data) throws NumberFormatException, IOException {
		Collection<String> merged=null;
		for(SingleInterval region: cluster.getAllIntervals()){
			Collection<String> barcodes=data.getBarcodesOverlappingRegion(region);
			merged=merge(merged, barcodes);
		}
	
		return merged;
	}

	private static Collection<String> merge(Collection<String> merged, Collection<String> barcodes) {
		if(merged==null){
			merged=new TreeSet<String>();
			merged.addAll(barcodes);
		}
		else{
			Collection<String> temp=new TreeSet<String>();
			for(String barcode: merged){
				if(barcodes.contains(barcode)){temp.add(barcode);}
			}
			merged=temp;
		}
		return merged;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		
		Cluster c=getClusterFromString(args[1]);
		System.err.println(c.toKmerString());
		
		FileWriter writer=new FileWriter(args[2]);
		
		int counter=0;
		while(data.hasNext()){
			Cluster cluster=data.next();
			
			if(cluster.getClusterSize()>1 && cluster.getClusterSize()<1000){
				if(cluster.containsAllIntervals(c.getAllIntervals())){
					writer.write(cluster.toSPRITEFormat()+"\n");
				}
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		data.close();
		writer.close();
	}

	private static Cluster getClusterFromString(String string) {
		Cluster rtrn=new Cluster("c");
		//chr10-0 chr22-19 chr10-132 :11
		//chr16-34475000_34500000 chr16-34525000_34550000 chr19-27725000_27750000 chr16-34575000_34600000 chr1-121475000_121500000 chr16-34650000_34675000 chr16-34700000_34725000
		
		String[] tokens=string.split(" ");
		
		for(int i=0; i<tokens.length; i++){
			String chr=tokens[i].split("-")[0];
			String region=tokens[i].split("-")[1];
			int start=new Integer(region.split("_")[0]);
			int end=new Integer(region.split("_")[1]);
			SingleInterval newInterval=new SingleInterval(chr, start, end);
			rtrn.addRead(newInterval);
		}
		return rtrn;
	}

	private static void write(String save, Collection<String> barcodes, BarcodingData data) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Cluster> clusters=data.getClustersWithBarcodes(barcodes);
		
		for(Cluster c: clusters){
			writer.write(c.toSPRITEFormat()+"\n");
		}
				
		writer.close();
	}
	
}
