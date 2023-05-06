package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class RNARNAMAp {

	
	//TODO Weight by number of ribosomes in cluster
	public RNARNAMAp(BarcodingDataStreaming data, Map<String, Kmer> rnaAliases, Kmer kmer, String save) throws IOException, InterruptedException{
		
		Map<String, Integer> rnaCounts=new TreeMap<String, Integer>();
		Map<String, Integer> totalCounts=new TreeMap<String, Integer>();
		Map<String, Integer> noDNACounts=new TreeMap<String, Integer>();
		Map<String, Integer> twoToTen=new TreeMap<String, Integer>();
		Map<String, Integer> tenToHundred=new TreeMap<String, Integer>();
		Map<String, Integer> hundredToThousand=new TreeMap<String, Integer>();
		Map<String, Double> weighted=new TreeMap<String, Double>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			if(c.containsRNA(kmer, true, rnaAliases)){
				increment(rnaCounts, c.getAllRNARegions());
				if(!c.hasDNA()){
					increment(noDNACounts, c.getAllRNARegions());
					weight(weighted, c, kmer, rnaAliases);
					if(c.getClusterSize()>1 && c.getClusterSize()<=10){increment(twoToTen, c.getAllRNARegions());}
					if(c.getClusterSize()>10 && c.getClusterSize()<=100){increment(tenToHundred, c.getAllRNARegions());}
					if(c.getClusterSize()>100 && c.getClusterSize()<1000){increment(hundredToThousand, c.getAllRNARegions());}
				}
			}
			increment(totalCounts, c.getAllRNARegions());
			counter++;
			if(counter%100000==0){System.err.println(counter);}
		}
		
		write(save, rnaCounts, totalCounts, noDNACounts, twoToTen, tenToHundred, hundredToThousand, weighted);
	}
	
	/*private void weight(Map<String, Double> weighted, Cluster c, Kmer kmer) {
		Collection<RNAInterval> rnas=c.getAllRNARegions();
		
		int num=0;
		for(RNAInterval rna: rnas){
			if(kmer.getRegions().contains(rna.getName())){num++;}
		}
		
		double weight=(double)num/(double)c.getRNAClusterSize();
		
		//System.err.println(num +"\t"+c.getClusterSize()+"\t"+c.getRNAClusterSize()+" "+weight);
		
		for(RNAInterval rna: rnas){
			double count=0;
			if(weighted.containsKey(rna.getName())){count=weighted.get(rna.getName());}
			count+=weight;
			weighted.put(rna.getName(), count);
		}
		
	}*/
	
	private void weight(Map<String, Double> weighted, Cluster c, Kmer kmer, Map<String, Kmer> aliases) {
		Collection<RNAInterval> rnas=c.getAllRNARegions();
		
		Kmer expanded=new Kmer();
		for(String rna: kmer.getRegions()){
			expanded.addRegions(get3(aliases,rna).getRegions());
		}
		
		int num=0;
		for(RNAInterval rna: rnas){
			if(expanded.getRegions().contains(rna.getName())){num++;}
		}
		
		double weight=(double)num/(double)c.getRNAClusterSize();
		
		//System.err.println(num +"\t"+c.getClusterSize()+"\t"+c.getRNAClusterSize()+" "+weight);
		
		for(RNAInterval rna: rnas){
			double count=0;
			if(weighted.containsKey(rna.getName())){count=weighted.get(rna.getName());}
			count+=weight;
			weighted.put(rna.getName(), count);
		}
		
	}

	private Kmer get3(Map<String, Kmer> collapseSets, String region) {
		if(collapseSets.containsKey(region)){return collapseSets.get(region);}
		
		Kmer rtrn=new Kmer();
		rtrn.addRegion(region);
		return rtrn;
	}

	private void weight(Map<String, Double> weighted, Cluster c) {
		double weight=1.0/(double)c.getClusterSize();
		
		Collection<RNAInterval> rnas=c.getAllRNARegions();
		for(RNAInterval rna: rnas){
			double count=0;
			if(weighted.containsKey(rna.getName())){count=weighted.get(rna.getName());}
			count+=weight;
			weighted.put(rna.getName(), count);
		}
		
	}

	private void increment(Map<String, Integer> rnaCounts, Collection<RNAInterval> rnas) {
		for(RNAInterval rna: rnas){
			int count=0;
			if(rnaCounts.containsKey(rna.getName())){count=rnaCounts.get(rna.getName());}
			count++;
			rnaCounts.put(rna.getName(), count);
		}
		
	}

	private void write(String save, Map<String, Integer> rnaCounts, Map<String, Integer> totalCounts, Map<String, Integer> noDNACounts, Map<String, Integer> twoToTen, Map<String, Integer> tenToHundred, Map<String, Integer> hundredToThousand, Map<String, Double> weighted) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("Gene\t Paired count \t Paired count without DNA \t All RNA Count\t 2-10 \t 10-100 \t 100-1000 \t weighted \n");
		
		for(String gene: totalCounts.keySet()){
			int totalcount=totalCounts.get(gene);
			int rnacount=get(rnaCounts, gene);
			int dnacount=get(noDNACounts, gene);
			int two=get(twoToTen, gene);
			int ten=get(tenToHundred, gene);
			int hundred=get(hundredToThousand, gene);
			double weight=get2(weighted, gene);
			writer.write(gene+"\t"+rnacount+"\t"+dnacount+"\t"+totalcount+"\t"+two+"\t"+ten+"\t"+hundred+"\t"+weight+"\n");
		}
		
		writer.close();
	}
	
	
	private double get2(Map<String, Double> weighted, String gene) {
		double rtrn=0;
		if(weighted.containsKey(gene)){rtrn=weighted.get(gene);}
		return rtrn;
	}

	private int get(Map<String, Integer> rnaCounts, String gene) {
		int rtrn=0;
		if(rnaCounts.containsKey(gene)){rtrn=rnaCounts.get(gene);}
		return rtrn;
	}
	
	private static Map<String, Kmer> parse(String string) throws IOException {
		Map<String, Kmer> rtrn=new TreeMap<String, Kmer>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String rna=tokens[0];
			String alias=tokens[1];
			Kmer kmer=new Kmer();
			if(rtrn.containsKey(alias)){kmer=rtrn.get(alias);}
			kmer.addRegion(rna);
			rtrn.put(alias, kmer);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save = args[1];
		Map<String, Kmer> aliases=parse(args[2]);
		
		Collection<String> names=new TreeSet<String>();
		for(int i=3; i<args.length; i++){
			names.add(args[i]);
		}
		
		Kmer kmer=new Kmer();
		kmer.addRegions(names);
		new RNARNAMAp(data, aliases, kmer, save);
		
	}else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=cluster file \n args[1]=output file name \n args[2]=alias file \n args[3-n]=RNA names";
	
}
