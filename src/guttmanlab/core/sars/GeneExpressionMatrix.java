package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GeneExpressionMatrix {

	
	public GeneExpressionMatrix(File bam, String save) throws IOException {
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<String, Integer> geneCount=new TreeMap<String, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String gene=read.getReferenceName();
				
				int score=0;
				if(geneCount.containsKey(gene)){score=geneCount.get(gene);}
				score++;
				geneCount.put(gene, score);
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		write(save, geneCount);
		
	}

	private void write(String save, TreeMap<String, Integer> geneCount) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: geneCount.keySet()) {
			writer.write(gene+"\t"+geneCount.get(gene)+"\n");
		}
		
		writer.close();
	}
	
}
