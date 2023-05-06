package guttmanlab.core.spidr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class NormalizedCounts {
	
	double total;
	
	public NormalizedCounts(File bam, int binSize, String save) throws IOException {
		Map<SingleInterval, Integer> scores=getScores(bam, binSize);
		write(scores, total, save);
	}
	
	public NormalizedCounts(File bam, int binSize, double freq, String save) throws IOException {
		Map<SingleInterval, Integer> scores=getScores(bam, binSize, freq);
		write(scores, total, save);
	}
	
	public NormalizedCounts(File[] bams, int binSize, List<String> regions, Map<String, Double> freqByName, String save) throws IOException {
		List<String> columns=getColumns(bams, freqByName);
		MatrixWithHeaders matrix=new MatrixWithHeaders(regions, columns);
		
		for(int i=0; i<bams.length; i++) {
			double freq=getFreq(freqByName, bams[i]);
			System.err.println(bams[i].getAbsolutePath()+"\t"+freq+"\t"+i+"\t"+bams.length);
			if(freq>0) {
				Map<SingleInterval, Integer> scores=getScores(bams[i], binSize);
				for(SingleInterval region: scores.keySet()) {
					int score=scores.get(region);
					if(matrix.containsRow(region.toUCSCStrand())) {
						matrix.set(region.toUCSCStrand(), bams[i].getName(), score);
					}
				}
			}
		}
		
		
		matrix.write(save);
	}

	
	private List<String> getColumns(File[] bams, Map<String, Double> freqByName) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<bams.length; i++) {
			double freq=getFreq(freqByName, bams[i]);
			if(freq>0) {rtrn.add(bams[i].getName());}
		}
		
		return rtrn;
	}

	private double getFreq(Map<String, Double> freqByName, File file) {
		if(freqByName.containsKey(file.getName())) {
			return freqByName.get(file.getName());
		}
		return -1;
	}

	private void write(Map<SingleInterval, Integer> scores, double total2, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		//writer.write("total\t"+total2+"\n");
		for(SingleInterval region: scores.keySet()) {
			int score=scores.get(region);
			//double norm=1000000*((double)score/total2);
			writer.write(region.toUCSCStrand()+"\t"+score+"\n");
		}
		writer.close();
	}

	private Map<SingleInterval, Integer> getScores(File bam, int binSize) {
		return getScores(bam, binSize, 1);
	}

	private Map<SingleInterval, Integer> getScores(File bam, int binSize, double freq) {
		SAMFileReader reader=new SAMFileReader(bam);
		
		SAMRecordIterator reads=reader.iterator();
		
		Map<SingleInterval, Integer> scores=new TreeMap<SingleInterval, Integer>();
		
		int counter=0;
		while(reads.hasNext()){
			double rand=Math.random();
			SAMRecord record=reads.next();
			if(rand<freq) {
				SAMFragment f=new SAMFragment(record);
				Collection<SingleInterval> allBins=f.getSingleInterval().allBins(binSize);
				for(SingleInterval binned: allBins) {
					binned.setOrientation(f.getOrientation());
					int score=0;
					if(scores.containsKey(binned)) {
						score=scores.get(binned);
					}
					score++;
					scores.put(binned, score);
				}
			}
				
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		total=counter;
		
		reader.close();
		reads.close();
		return scores;
	}
	
	private static Map<String, Double> parse(String string) throws IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines) {
			rtrn.put(line.split("\t")[0], Double.parseDouble(line.split("\t")[2]));
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			/*File[] bams=new File(args[0]).listFiles();
			int binSize=100;
			Map<String, Double> freqByName=parse(args[1]);
			List<String> regions=BEDFileIO.loadLines(args[2]);
			String save=args[3];
			new NormalizedCounts(bams, binSize, regions, freqByName, save);*/
			
			
			File bam=new File(args[0]);
			String save=args[1];
			double freq=Double.parseDouble(args[2]);
			int binSize=100;
			new NormalizedCounts(bam, binSize, freq, save);
			
		}else {System.err.println(usage);}
		
	}

	static String usage=" args[0]=bam \n args[1]=save \n args[2]=sample frequency";
	
}
