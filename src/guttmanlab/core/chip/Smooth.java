package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class Smooth {

	
	Map<SingleInterval, Integer> startCounts;
	Map<SingleInterval, Integer> endCounts;
	int binSize;
	
	
	public Smooth(File file, int binSize, String save) throws IOException {
		this.binSize=binSize;
		Map<SingleInterval, Integer> counts=scoreWindows(file, binSize);
		
		
		//binSize/2 up and down
		/*Map<SingleInterval, Integer> map=new TreeMap<SingleInterval, Integer>();
		int counter=0;
		for(SingleInterval r: counts.keySet()) {
			Map<SingleInterval, Integer> list1=extendUp(r, counts.get(r), binSize); // add new starts, subtract old ends
			Map<SingleInterval, Integer> list2=extendDown(r, counts.get(r), binSize); //add new ends, subtract old starts
			map.putAll(list1);
			map.putAll(list2);
			counter++;
			if(counter%10000==0) {System.err.println(counter+" "+counts.size());}
		}*/
		
		
		counts=filter(counts, "chr7");
		write(save, counts);
		
	}
	
	private void write(String save, Map<SingleInterval, Integer> counts) throws IOException {
		//Map<String, Integer> chrSizes=CoordinateSpace.getGenomeLengths(fileHeader);
		FileWriter writer=new FileWriter(save);
		
		
		//Go through starts and bin around
		int counter=0;
		
		for(SingleInterval r: counts.keySet()) {
			int start=r.getReferenceStartPosition();
			int end=r.getReferenceEndPosition();
			for(int i=start; i<end; i++) {
				int newStart=start-binSize/2;
				double sum=getScore(r.getReferenceName(), newStart, newStart+binSize);
				writer.write(r.getReferenceName()+"\t"+start+"\t"+start+"\t"+sum+"\n");
			}
			counter++;
			if(counter%10000==0) {System.err.println(counter+" "+this.startCounts.size());}
		}
		
		/*for(SingleInterval r: this.startCounts.keySet()) {
			double sum=getScore(r.getReferenceName(), r.getReferenceStartPosition()-binSize/2, r.getReferenceEndPosition()+binSize/2);
			writer.write(r.toBedgraph(sum)+"\n");
			counter++;
			if(counter%10000==0) {System.err.println(counter+" "+this.startCounts.size());}
		}*/
		
		/*for(String chr: chrSizes.keySet()) {
			System.err.println(chr);
			for(int i=0; i<chrSizes.get(chr); i++) {
				double score=getScore(chr, i);
				if(score>0) {writer.write(chr+"\t"+i+"\t"+(i+1)+"\t"+score+"\n");}
			}
		}*/
		writer.close();
	}

	private double getScore(String chr, int start, int end) {
		double sum=0;
		for(int i=start; i<end; i++) {
			sum+=get(this.startCounts, new SingleInterval(chr, i, i));
		}
		return sum;
	}

	private Map<SingleInterval, Integer> filter(Map<SingleInterval, Integer> counts, String chr) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		for(SingleInterval r: counts.keySet()) {
			if(r.getReferenceName().equals(chr)) {
				rtrn.put(r, counts.get(r));
			}
		}
		
		return rtrn;
	}

	/*private void write(String save, Map<SingleInterval, Integer> map) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: map.keySet()) {
			int score=map.get(r);
			//writer.write(r.getReferenceName()+"\t"+r.getReferenceStartPosition()+"\t"+(r.getReferenceStartPosition()+1)+"\t"+score+"\n");
			writer.write(r.getMidPoint().toBedgraph(score)+"\n");
		}
		
		writer.close();
	}*/

	private Map<SingleInterval, Integer> extendUp(SingleInterval r, int baseScore, int binSize) throws IOException {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		for(int i=0; i<binSize/2; i++) {
			SingleInterval newInterval= new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition()-i, r.getReferenceEndPosition()-i);
			//System.err.println(i+" "+r.toUCSC()+" "+newInterval.toUCSC());
			//add starts
			int add=startCounts(r.getReferenceName(), newInterval.getReferenceStartPosition(), r.getReferenceStartPosition());
			int subtract=endCounts(r.getReferenceName(), newInterval.getReferenceEndPosition(), r.getReferenceEndPosition());
			int newScore=baseScore+add-subtract;
			rtrn.put(newInterval, newScore);
		}
		return rtrn;
	}
	
	
	private Map<SingleInterval, Integer> extendDown(SingleInterval r, int baseScore, int binSize) throws IOException {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		for(int i=0; i<binSize/2; i++) {
			SingleInterval newInterval= new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition()+i, r.getReferenceEndPosition()+i);
			
			//add starts
			int subtract=startCounts(r.getReferenceName(), r.getReferenceStartPosition(), newInterval.getReferenceStartPosition());
			int add=endCounts(r.getReferenceName(), r.getReferenceEndPosition(), newInterval.getReferenceEndPosition());
			int newScore=baseScore+add-subtract;
			rtrn.put(newInterval, newScore);
		}
		return rtrn;
	}

	private int endCounts(String chr, int start, int end) {
		int sum=0;
		for(int i=start; i<end; i++) {
			SingleInterval pos=new SingleInterval(chr, i, i);
			sum+=get(endCounts,pos);
		}
		return sum;
	}

	private int startCounts(String chr, int start, int end) {
		int sum=0;
		for(int i=start; i<end; i++) {
			SingleInterval pos=new SingleInterval(chr, i, i);
			sum+=get(startCounts,pos);
		}
		return sum;
	}

	private int get(Map<SingleInterval, Integer> startCounts2, SingleInterval pos) {
		int rtrn=0;
		if(startCounts2.containsKey(pos)) {rtrn=startCounts2.get(pos);}
		return rtrn;
	}

	private Map<SingleInterval, Integer> scoreWindows(File file, int binSize) {
		this.startCounts=new TreeMap<SingleInterval, Integer>();
		this.endCounts=new TreeMap<SingleInterval, Integer>();
		
		
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			addEnds(record);
			SingleInterval binned=SingleInterval.bin(record, binSize);
			int score=0;
			if(counts.containsKey(binned)) {score=counts.get(binned);}
			score++;
			counts.put(binned, score);
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reader.close();
		reads.close();
		return counts;	
	}
	
	private void addEnds(SAMRecord record) {
		addStart(record);
		addEnd(record);
		
	}


	private void addStart(SAMRecord record) {
		SingleInterval r=new SingleInterval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentStart());
		int count=0;
		if(this.startCounts.containsKey(r)) {count=this.startCounts.get(r);}
		count++;
		this.startCounts.put(r, count);
		
	}


	private void addEnd(SAMRecord record) {
		SingleInterval r=new SingleInterval(record.getReferenceName(), record.getAlignmentEnd(), record.getAlignmentEnd());
		int count=0;
		if(this.endCounts.containsKey(r)) {count=this.endCounts.get(r);}
		count++;
		this.endCounts.put(r, count);
	}
	
	public static void main (String[] args) throws IOException {
		File bam1=new File(args[0]);
			String save=args[1];
			int binSize=Integer.parseInt(args[2]);
			new Smooth(bam1, binSize, save);
	}
	
}
