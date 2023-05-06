package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SingleInterval;

public class ExonIntronBedgraph {
	
	int resolution=1000000;
	
	static SingleInterval problemBin1=new SingleInterval("chr2:79490000-79500000");
	static SingleInterval problemBin2=new SingleInterval("chr11:3119270-3192250");
	static SingleInterval problemBin3=new SingleInterval("chr15:99734977-99736026");
	static SingleInterval problemBin4=new SingleInterval("chr3:5173978-5175025");
	static SingleInterval problemBin5=new SingleInterval("chr13:58176952-58178051");
	
	static Collection<SingleInterval> problemBins=new TreeSet<SingleInterval>();


	public ExonIntronBedgraph(BarcodingDataStreaming data, String gene, String type, String save, int binResolution, boolean weight) throws IOException{
		
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		Collection<Cluster> clusters=getRNAClusters(data, gene, type);
		
		double total=0;
		double onDNA=0;
		for(Cluster c: clusters){
			if(c.getAllDNAIntervals().size()>0){onDNA++;}
			total++;
			
			Cluster binned=c.bin(binResolution);
				for(SingleInterval region: binned.getAllDNAIntervals()){
					double score=0;
					if(rtrn.containsKey(region)){
						score=rtrn.get(region);
					}
				if(weight){score+=(2.0/binned.getAllDNAIntervals().size());}
				else{score++;}
				rtrn.put(region, score);
			}
		}
			
		double fraction=onDNA/total;
		System.err.println(onDNA+" "+total+" "+fraction);
		
		
		rtrn=removeProblemBins(rtrn);
		data.close();
		write(save, rtrn, total-onDNA);
	}
	
	
	private Collection<Cluster> getRNAClusters(BarcodingDataStreaming data, String gene, String type) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		while(data.hasNext()){
			Cluster c=data.next();
			boolean useCluster=useCluster(c, gene, type);
			if(useCluster){rtrn.add(c);}
		}
		
		data.close();
		
		return rtrn;
	}


	private boolean useCluster(Cluster c, String gene, String type) {
		for(RNAInterval rna: c.getAllRNARegions()){
			Collection<String> rnas=rna.getExonOnlyNames();
			if(type.equalsIgnoreCase("intron")){rnas=rna.getIntronOnlyNames();}
			else if(type.equalsIgnoreCase("all") || type.equalsIgnoreCase("both")){
				rnas=rna.getAllRNANames();
			}
			
			if(rnas.contains(gene)){return true;}
		}
		return false;
	}


	private void write(String save, Map<SingleInterval, Double> map, double offDNA) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: map.keySet()){
			double val=map.get(region);
			val=100*(val/offDNA);
			writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+val+"\n");
		}
		
		writer.close();
	}
	
	
	private static boolean overlaps(SingleInterval region, Collection<SingleInterval> problemBins) {
		for(SingleInterval problemBin: problemBins){
			if(region.overlaps(problemBin)){return true;}
		}
		return false;
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

	
	
	private boolean contains(Collection<RNAInterval> rnas, String gene, String type) {
		for(RNAInterval rna: rnas){
			if(rna.getName().equals(gene)){
				if(type.equalsIgnoreCase("Exon")){return rna.isExon();}
				else if(type.equalsIgnoreCase("Intron")){return rna.isIntron();}
			}
		}
		return false;
	}


	
	
	public static void main(String[] args) throws IOException{
		problemBins.add(problemBin1);
		problemBins.add(problemBin2);
		problemBins.add(problemBin3);
		problemBins.add(problemBin4);
		problemBins.add(problemBin5);
		
		if(args.length>3){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String gene=args[1];
		String type=args[2];
		String save=args[3];
		
		int binResolution=1000000;
		boolean weight=false;
		new ExonIntronBedgraph(data, gene, type, save, binResolution, weight);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]= data \n args[1]=gene \n args[2]=type \n args[3]=save";
}
