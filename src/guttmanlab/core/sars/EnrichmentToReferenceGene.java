package guttmanlab.core.sars;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class EnrichmentToReferenceGene {

	public EnrichmentToReferenceGene(File file1, File file2, String testGene, String referenceGene) throws IOException{
		TreeMap<String, Double> sample=score(file1);
		TreeMap<String, Double> input=score(file2);
		
		double fractionCLAP=sample.get(testGene)/input.get(testGene);
				
		double fractionInput=sample.get(referenceGene)/input.get(referenceGene);
		double ratio=fractionCLAP/fractionInput;
		System.err.println(fractionCLAP+" "+fractionInput+" "+ratio);
		
		
	}
	
	
	
	
	private TreeMap<String, Double> score(File bam){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<String, Double> positionCount=new TreeMap<String, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String chr=read.getReferenceName();
				double score=0;
				if(positionCount.containsKey(chr)){score=positionCount.get(chr);}
				score++;
				positionCount.put(chr, score);
				totalCount++;
			}
				
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		
		
		return positionCount;
	}
	
	public static void main(String[] args) throws IOException{
		File file1=new File(args[0]);
		File file2=new File(args[1]);
		String testGene=args[2];
		String referenceGene=args[3];
		new EnrichmentToReferenceGene(file1, file2, testGene, referenceGene);
	}
	
}
