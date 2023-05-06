package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;

public class FractionSplicedFromSPRITE {
	
	int binResolution=1000000;

	public FractionSplicedFromSPRITE(BarcodingDataStreaming data, String save) throws IOException{
		
		FileWriter writer =new FileWriter(save);
		
		//Get RNA on it's DNA locus
		getNascent(data, writer);
		
		
		writer.close();

		
		//Count spliced versus unspliced
		
	}
	
	
	private void getNascent(BarcodingDataStreaming data, FileWriter writer) throws IOException {
		int counter=0;
		int total=0;
		while(data.hasNext()){
			Cluster c=data.next();
			int size=getNascent(c, writer);
			total+=size;
			counter++;
			if(counter%100000 ==0){System.err.println(counter+" "+total);}
		}
		
		data.close();
	}

	private int getNascent(Cluster c, FileWriter writer) throws IOException {
		Collection<RNAInterval> rtrn=new TreeSet<RNAInterval>();
		
		Cluster binned=c.bin(binResolution);
		
		for(RNAInterval rnaRegion:c.getAllRNARegions()){
			SingleInterval rnaBin=rnaRegion.bin(binResolution);
			for(SingleInterval dnaRegion: binned.getAllDNAIntervals()){
				if(rnaBin.overlaps(dnaRegion)){rtrn.add(rnaRegion);}
			}
		}
		
		int size=0;
		for(RNAInterval region: rtrn){
			if(region.getType().equalsIgnoreCase("exon") || region.getType().equalsIgnoreCase("intron")){
				writer.write(region.toBED()+"\t"+region.getType()+"\n");
				size++;
			}
		}
		return size;
	}

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		new FractionSplicedFromSPRITE(data, save);
	}
	
}
