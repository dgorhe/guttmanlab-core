package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;

public class MakeRNADNAContactBedGraphMultiple {

	/*static SingleInterval problemBin1=new SingleInterval("chr2:79490000-79500000");
	static SingleInterval problemBin2=new SingleInterval("chr11:3119270-3192250");
	static SingleInterval problemBin3=new SingleInterval("chr15:99734977-99736026");
	static SingleInterval problemBin4=new SingleInterval("chr3:5173978-5175025");
	static SingleInterval problemBin5=new SingleInterval("chr13:58176952-58178051");*/
	
	static Collection<SingleInterval> problemBins=new TreeSet<SingleInterval>();

	
	private static void write(String save, Map<SingleInterval, Double> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			double val=map.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		
		writer.close();
	}
	
	private static void write(String save, Map<SingleInterval, Double> map, double minCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		double max=max(map);
		
		for(SingleInterval region: map.keySet()){
			double val=map.get(region);
			//if(val>=minCount) {
			double score=100*(val/max);
				writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+score+"\n");
			//}
		}
		
		writer.close();
	}
	
	
	
	private static double max(Map<SingleInterval, Double> map) {
		return Statistics.max(map.values());
	}

	private static void writeSmoothed(String save, Map<SingleInterval, Double> map, int binResolution, int n) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			Collection<SingleInterval> regions=getNeighbors(region, binResolution, n);
			double val=average(regions, map);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		writer.close();
	}
	
	private static Collection<SingleInterval> getNeighbors(SingleInterval region, int binResolution, int n) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		//upstream
		for(int i=0; i<n; i++){
			int start=region.getReferenceStartPosition()-((i+1)*binResolution);
			int end=start+binResolution;
			SingleInterval newRegion=new SingleInterval(region.getReferenceName(), start, end);
			rtrn.add(newRegion);
		}
		
		
		//downstream
		for(int i=0; i<n; i++){
			int start=region.getReferenceStartPosition()+((i+1)*binResolution);
			int end=start+binResolution;
			SingleInterval newRegion=new SingleInterval(region.getReferenceName(), start, end);
			rtrn.add(newRegion);
		}
		
		return rtrn;
	}



	private static double average(Collection<SingleInterval> regions, Map<SingleInterval, Double> map) {
		double sum=0;
		double count=0;
		
		for(SingleInterval region: regions){
			if(map.containsKey(region)){
				sum+=map.get(region);
				count++;
			}
		}
		
		return sum/count;
	}

	private static Map<SingleInterval, Double> removeProblemBins(Map<SingleInterval, Double> map) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval region: map.keySet()){
			if(!overlaps(region, problemBins)){
				rtrn.put(region, map.get(region));
			}
		}
		
		return rtrn;
	}


	private static boolean overlaps(SingleInterval region, Collection<SingleInterval> problemBins) {
		for(SingleInterval problemBin: problemBins){
			if(region.overlaps(problemBin)){return true;}
		}
		return false;
	}

	
	private static Map<SingleInterval, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		Map<SingleInterval, String> rtrn=new TreeMap<SingleInterval, String>();
		for(String line: lines){
			rtrn.put(new SingleInterval(line.split("\t")[0]), line.split("\t")[1]);
		}
		return rtrn;
	}

	
	private static Collection<String> getClasses(Map<SingleInterval, String> genes) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(SingleInterval rna: genes.keySet()){
			String rnaClass=genes.get(rna);
			rtrn.add(rnaClass);
		}
		
		return rtrn;
	}
	
	private static Kmer getKmer(Map<SingleInterval, String> genes, String rnaClass) {
		Kmer rtrn=new Kmer();
		for(SingleInterval rna: genes.keySet()){
			String rnaClass2=genes.get(rna);
			if(rnaClass2.equals(rnaClass)){rtrn.addInterval(rna);}
		}
		return rtrn;
	}
	
	private static Map<SingleInterval, Double> getCounts(Collection<SingleInterval> regions, int binResolution) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval r: regions) {
			Collection<SingleInterval> bins=r.allBins(binResolution);
			for(SingleInterval bin: bins) {
				double count=0;
				if(rtrn.containsKey(bin)) {count=rtrn.get(bin);}
				count++;
				rtrn.put(bin, count);
			}
		}
		
		return rtrn;
	}
	
	private static void write(String save, Map<SingleInterval, Double> counts, Map<SingleInterval, Double> clusterCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		boolean isEmpty=true;
		
		for(SingleInterval region: counts.keySet()) {
			double readCount=counts.get(region);
			double clusterCount=0;
			if(clusterCounts.containsKey(region)) {clusterCount=clusterCounts.get(region);}
			if(readCount>1 && clusterCount>1) {
				double min=Math.min(readCount, clusterCount);
				writer.write(region.toBedgraph(min)+"\n");
				isEmpty=false;
				//writer.write(region.toBedgraph(clusterCount)+"\n");
			}
		}
		
		writer.close();
		
		if(isEmpty) {
			System.err.println(save+" deleted");
			new File(save).delete();
		}
	}


	public static void main(String[] args) throws IOException{
		if(args.length>3){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Collection<String> geneNames=BEDFileIO.loadLines(args[1]);
			String saveDir=args[2];
			int binResolution=Integer.parseInt(args[3]);
			
			//data.getClusters(geneNames, saveDir);
			
			for(String gene: geneNames) {
				System.out.println(gene);
				String save=saveDir+"/"+gene+".counts.bedgraph";
				Map<SingleInterval, Double> counts=getCounts(data.getRegions(gene), binResolution);
				Map<SingleInterval, Double> clusterCounts=data.getRNADNAClusterCounts(binResolution, gene);
				
				write(save, counts, clusterCounts);
				
			}
		
			data.close();
			
		}
		else{System.err.println(usage);}
	}

	

	


	






	static String usage=" args[0]=clusters \n args[1]=gene list \n args[2]=save \n args[3]=binResolution";
	
}
