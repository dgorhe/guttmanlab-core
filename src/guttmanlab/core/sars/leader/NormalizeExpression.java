package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class NormalizeExpression {

	public NormalizeExpression(File rRNAFile, File ncFile, File mRNAFile, String save) throws IOException {
		Map<String, Double> rRNAScores=getScores(rRNAFile);
		Map<String, Double> ncRNAScores=getScores(ncFile);
		Map<String, Double> mRNAScores=getScores(mRNAFile);
		
		FileWriter writer=new FileWriter(save);
		double avg=average(rRNAScores);
		
		for(String gene: ncRNAScores.keySet()) {
			double score=ncRNAScores.get(gene);
			double ratio=score/avg;
			writer.write(gene+"\t"+"ncRNA\t"+score+"\t"+ratio+"\n");
		}
		
		
		for(String gene: mRNAScores.keySet()) {
			double score=mRNAScores.get(gene);
			double ratio=score/avg;
			writer.write(gene+"\t"+"mRNA\t"+score+"\t"+ratio+"\n");
		}
		
		writer.close();
	}
	
	
	private double average(Map<String, Double> rRNAScores) {
		double sum=0;
		
		for(String rna: rRNAScores.keySet()) {sum+=rRNAScores.get(rna);}
		
		return sum/(double)rRNAScores.size();
	}


	private Map<String, Double> getScores(File file) {
		SAMFileReader inputReader= new SAMFileReader(file);
		Map<String, Double> geneCount=new TreeMap<String, Double>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				String gene=read.getReferenceName();
				
				double score=0;
				if(geneCount.containsKey(gene)){score=geneCount.get(gene);}
				score++;
				geneCount.put(gene, score);
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		return geneCount;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			File file1=new File(args[0]);
			File file2=new File(args[1]);
			File file3=new File(args[2]);
			String save=args[3];
			new NormalizeExpression(file1, file2, file3, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=45S \n args[1]=noncoding \t args[2]=mRNA \n args[3]=save";
	
}
