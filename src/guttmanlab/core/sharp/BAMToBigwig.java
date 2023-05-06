package guttmanlab.core.sharp;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BAMToBigwig {

	public BAMToBigwig(File bam, String save, int binSize) throws IOException, InterruptedException {
		Map<SingleInterval, Double> scores=score(bam, binSize);
		BEDFileIO.writeBEDGraph(scores, save);
		//toBigwig(save, sizes);
	}
	
	
	private void toBigwig(String save, String sizesFile) throws IOException, InterruptedException {
		String cmd="/central/groups/guttman/software/kentUtils/bin/bedGraphToBigWig "+save+" "+sizesFile+" "+save+".bw";
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
	}


	private Map<SingleInterval, Double> score(File bam, int binSize){
		SAMFileReader inputReader= new SAMFileReader(bam);
		TreeMap<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && read.getMappingQuality()>1){
				Collection<SingleInterval> allBins=SAMFragment.allBins(read, binSize);
				for(SingleInterval bin: allBins) {
					//System.out.println(bin);
					int score=0;
					if(positionCount.containsKey(bin)){score=positionCount.get(bin);}
					score++;
					positionCount.put(bin, score);
				}
			}
				totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
				
		reads.close();
		inputReader.close();
		
		
		Map<SingleInterval, Double> norm=new TreeMap<SingleInterval, Double>();
		
		for(SingleInterval r: positionCount.keySet()) {
			int score=positionCount.get(r);
			double val=Math.pow(10, 7) *((double)score/(double)totalCount);
			norm.put(r, val);
		}
		
		return norm;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		File bam=new File(args[0]);
		String save=args[1];
		//String sizes=args[2];
		new BAMToBigwig(bam, save, 10);
	}
	
}
