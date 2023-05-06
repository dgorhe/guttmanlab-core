package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class TranslationalDecay {

	
	
	public TranslationalDecay(File bamFile1,  String save) throws IOException{
		
		Map<String, Double> countsByGene1=score(bamFile1);
		
		write(save, countsByGene1);
		
	}
	
	public TranslationalDecay(File bamFile1, File bamFile2, String save) throws IOException{
		
		Map<String, Double> countsByGene1=score(bamFile1);
		Map<String, Double> countsByGene2=score(bamFile2);
		
		write(save, countsByGene1, countsByGene2);
		
	}

	private Map<String, Double> score(File bamFile1) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		SAMFileReader inputReader= new SAMFileReader(bamFile1);
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String gene=read.getReferenceName();
			double score=0;
			if(rtrn.containsKey(gene)){score=rtrn.get(gene);}
			score+=1.0;
			rtrn.put(gene, score);
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		Map<String, Double> newRtrn=new TreeMap<String, Double>();
		for(String gene: rtrn.keySet()){
			double val=rtrn.get(gene);
			double newVal=1000000*(val/(double)totalCount);
			newRtrn.put(gene, newVal);
		}
		
		return newRtrn;
	}

	private void write(String save, Map<String, Double> countsByGene1, Map<String, Double> countsByGene2) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: countsByGene1.keySet()){
			double val1=get(countsByGene1, gene);
			double val2=get(countsByGene2, gene);
			double ratio=val1/val2;
			writer.write(gene+"\t"+val1+"\t"+val2+"\t"+ratio+"\n");
		}
		
		writer.close();
	}
	
	private void write(String save, Map<String, Double> countsByGene1) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: countsByGene1.keySet()){
			double val1=get(countsByGene1, gene);
			writer.write(gene+"\t"+val1+"\n");
		}
		
		writer.close();
	}

	private double get(Map<String, Double> countsByGene1, String gene) {
		if(countsByGene1.containsKey(gene)){return countsByGene1.get(gene);}
		return 0;
	}
	
	public static void main(String[] args) throws IOException{
		File file1=new File(args[0]);
		String save=args[1];
		new TranslationalDecay(file1, save);
	}
	
}
