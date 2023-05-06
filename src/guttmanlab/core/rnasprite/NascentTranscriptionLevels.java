package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;

public class NascentTranscriptionLevels {

	
	public NascentTranscriptionLevels(BarcodingDataStreaming data, String save, int resolution) throws IOException{
		
		Map<SingleInterval, Integer> posScores=new TreeMap<SingleInterval, Integer>();
		//Map<SingleInterval, Integer> negScores=new TreeMap<SingleInterval, Integer>();
		
		
		
		//FileWriter writer= new FileWriter(save);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			Collection<RNAInterval> rnas=c.getAllRNARegions();
			for(RNAInterval rna: rnas){
				if(rna.getType().equals("intron")){
					String name=rna.getName();
					
					SingleInterval binned=rna.bin(resolution);
					add(posScores, binned);
					//else if(rna.getOrientation().equals(Strand.NEGATIVE)){add(negScores, binned);}
					
					//writer.write(rna.getReferenceName()+"\t"+rna.getReferenceStartPosition()+"\t"+rna.getReferenceEndPosition()+"\t"+name+"\n");
					
					/*int score=0;
					if(scores.containsKey(name)){score=scores.get(name);}
					score++;
					scores.put(name, score);*/
				}
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		write(save, posScores);
		//write(save+".negScores.bedgraph", negScores);
		
	}
	
	private void add(Map<SingleInterval, Integer> scores, SingleInterval name) {
		int score=0;
		if(scores.containsKey(name)){score=scores.get(name);}
		score++;
		scores.put(name, score);
	}

	private void write(String save, Map<SingleInterval, Integer> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: scores.keySet()){
			writer.write(region.tobedgraph(scores.get(region))+"\n");}
		
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		int resolution=new Integer(args[2]);
		new NascentTranscriptionLevels(data, save, resolution);
	}
	
}
