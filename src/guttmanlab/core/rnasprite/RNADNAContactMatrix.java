package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
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
import guttmanlab.core.datastructures.Pair;

public class RNADNAContactMatrix {

	static SingleInterval problemBin1=new SingleInterval("chr2:79490000-79500000");
	static SingleInterval problemBin2=new SingleInterval("chr11:3119270-3192250");
	static SingleInterval problemBin3=new SingleInterval("chr15:99734977-99736026");
	static SingleInterval problemBin4=new SingleInterval("chr3:5173978-5175025");
	static SingleInterval problemBin5=new SingleInterval("chr13:58176952-58178051");
	
	static Collection<SingleInterval> problemBins=new TreeSet<SingleInterval>();

	
	private static void write(String save, Map<SingleInterval, Double> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			double val=map.get(region);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		
		writer.close();
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
	
	private static List<String> removeProblemBins(List<String> list) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String r: list){
			SingleInterval region=new SingleInterval(r);
			if(!overlaps(region, problemBins)){
				rtrn.add(r);
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
	
	private static Map<String, Kmer> parse(String string) throws IOException {
		Map<String, Kmer> rtrn=new TreeMap<String, Kmer>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String group=tokens[1];
			String name=tokens[0].replaceAll("\"", "");
			Kmer kmer=new Kmer();
			if(rtrn.containsKey(group)){
				kmer=rtrn.get(group);
			}
			kmer.addRegion(name);
			rtrn.put(group, kmer);
		}
		
		return rtrn;
	}

	



	private static void add(String gene, Map<SingleInterval, Double> map, MatrixWithHeaders mwh) {
		for(SingleInterval region: map.keySet()){
			double score=map.get(region);
			mwh.set(gene, region.toUCSC(), score);
		}
		
	}



	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			Map<String, Kmer> genes=parse(args[1]);
			int binResolution=new Integer(args[2]);
			String save=args[3];
			boolean weightByClusterSize=new Boolean(args[4]);
			
			
			problemBins.add(problemBin1);
			problemBins.add(problemBin2);
			problemBins.add(problemBin3);
			problemBins.add(problemBin4);
			problemBins.add(problemBin5);
			
			
			int minCluster=0;
			int maxCluster=1000;
			
			List<String> columns=data.getGenomePositions(binResolution);
			columns=removeProblemBins(columns);
			
			List<String> rows=new ArrayList<String>(genes.keySet());
			
			MatrixWithHeaders mwh=new MatrixWithHeaders(rows, columns);
			for(String gene: genes.keySet()){
				Kmer kmer=genes.get(gene);
				System.err.println(gene +" "+kmer.getRegions());
				Map<SingleInterval, Double> map=data.getRNADNAContacts(binResolution, kmer, weightByClusterSize, false, minCluster, maxCluster);
				map=removeProblemBins(map);
				
				add(gene, map, mwh);
				
				data.close();
			}
			
			mwh.write(save);
		}
		else{System.err.println(usage);}
	}

	
	static String usage=" args[0]=clusters \n args[1]=gene class file \n args[2]=bin resolution \n args[3]=save \n args[4]=weight by cluster size (true or false)";
	
}
