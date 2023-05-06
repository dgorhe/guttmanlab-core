package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.annotation.Annotation.Strand;

public class SplicingRatios {
	
	int binResolution=1000;

	
	public SplicingRatios(BarcodingDataStreaming data, String save) throws IOException{
		
		//FileWriter writer=new FileWriter(save+".bed");
		//Find RNA that is paired with DNA in same region
		//Count exons versus introns
		
		//Map<String, Pair<Integer>> intronExonCounts=new TreeMap<String, Pair<Integer>>();
		
		Map<SingleInterval, Double> counts=new TreeMap<SingleInterval, Double>();
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Collection<RNAInterval> rnas=c.getAllRNARegions();
			for(RNAInterval rna: rnas){
				boolean isNascent=isNascent(rna, c);
				if(isNascent){
					Collection<SingleInterval> bins=rna.getBins(binResolution);
					for(SingleInterval bin: bins) {
						double count=0;
						if(counts.containsKey(bin)) {count=counts.get(bin);}
						count++;
						counts.put(bin, count);
					}
					//writer.write(rna.getReferenceName()+"\t"+rna.getReferenceStartPosition()+"\t"+rna.getReferenceEndPosition()+"\t"+rna.getType()+"\n");
					//add(rna, intronExonCounts);
				}
			}
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		
		BEDFileIO.writeBEDGraph(counts, save);
		
		//writer.close();
		//write(save, intronExonCounts);
	}
	
	private boolean isNascent(RNAInterval rna, Cluster c) {
		//check whether the cluster has a DNA overlapping the rna interval
		SingleInterval binned=rna.bin(binResolution);
		return c.containsOverlappingDNA(binned);
		
	}

	private void add(RNAInterval rna, Map<String, Pair<Integer>> intronExonCounts) {
		Pair<Integer> vals=new Pair<Integer>(0,0);
		
		if(intronExonCounts.containsKey(rna.getName())){vals=intronExonCounts.get(rna.getName());}
		
		int intronCount=vals.getValue1();
		int exonCount=vals.getValue2();
		if(rna.getType().equalsIgnoreCase("intron")){intronCount++;}
		if(rna.getType().equalsIgnoreCase("exon")){exonCount++;}
		
		vals.setValue1(intronCount);
		vals.setValue2(exonCount);
		
		intronExonCounts.put(rna.getName(), vals);
		
	}

	private void write(String save, Map<String, Pair<Integer>> intronExonCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String geneName: intronExonCounts.keySet()){
			Pair<Integer> vals=intronExonCounts.get(geneName);
			writer.write(geneName+"\t"+vals.getValue1()+"\t"+vals.getValue2()+"\n");
		}
		
		writer.close();
	}
	

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new SplicingRatios(data, save);
	}
	
}
