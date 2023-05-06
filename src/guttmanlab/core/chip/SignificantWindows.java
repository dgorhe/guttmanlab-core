package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.PoissonDistribution;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.math.MaximumContiguousSubsequence;
import guttmanlab.core.math.ScanStat;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SignificantWindows {

	double totalReads;
	double genomeLength;
	double numPositionsCovered;
	static double alpha=0.01;
	double threshold=0.05;
	double lambda;
	int binSize;
	Map<SingleInterval, Integer> startCounts;
	Map<SingleInterval, Integer> endCounts;
	int cutoff;
	
	
	public SignificantWindows(File file, int binSize, String save) throws IOException {
		this.binSize=binSize;
		
		Map<SingleInterval, Integer> sigRegions=significantWindows(file, binSize);
		
		Map<SingleInterval, Integer> collapsed=collapse(sigRegions);
		
		//Extend if needed
		Collection<SingleInterval> extend=extend(collapsed, sigRegions);
		
		//Trim ends
		Collection<SingleInterval> trimmed=trim(file, extend);
		
		write(save, trimmed);
	}
	
	private Collection<SingleInterval> extend(Map<SingleInterval, Integer> collapsed, Map<SingleInterval, Integer> sigRegions) {
		Collection<SingleInterval> extensions=new TreeSet<SingleInterval>();
		
		int counter=0;
		for(SingleInterval r: collapsed.keySet()) {
			String chr=r.getReferenceName();
			SingleInterval start=new SingleInterval(chr, r.getReferenceStartPosition(), r.getReferenceStartPosition()+binSize);
			SingleInterval end=new SingleInterval(chr, r.getReferenceEndPosition()-binSize, r.getReferenceEndPosition());
			int startScore=sigRegions.get(start);
			int endScore=sigRegions.get(end);
			
			int startExtension=getStartExtension(r, startScore);
			int endExtension=getEndExtension(r, endScore);
			
			SingleInterval newInterval=new SingleInterval(r.getReferenceName(), startExtension, endExtension);
			extensions.add(newInterval);
			counter++;
			if(counter%1000==0) {System.err.println(counter+" "+collapsed.size()+" "+cutoff);}
		}
		
		return collapse(extensions);
		//return extensions;
	}

	
	
	private int getEndExtension(SingleInterval r,  int score) {
		String chr=r.getReferenceName();
		int end=r.getReferenceEndPosition();
		int updatedScore=score;
		
		int maxEnd=r.getReferenceEndPosition();
		
		for(int i=1; i<=binSize; i++) {
			int newEnd=end+i;
			int newStart=newEnd-binSize;
			
			updatedScore=updatedScore-getCounts(chr, newStart-1, this.endCounts)+getCounts(chr, newEnd, this.startCounts);
			
			if(updatedScore>cutoff) {
				maxEnd=Math.max(maxEnd, newEnd);
			}
			
			/*int w=binSize;
			double pval=ScanStat.getPValue(updatedScore, lambda, w, genomeLength);
			if(pval<alpha) {
				maxEnd=Math.max(maxEnd, newEnd);
			}*/
		}
		
		return maxEnd;
		
		
		
		
		/*String chr=r.getReferenceName();
		int end=r.getReferenceEndPosition();
		for(int i=binSize; i>0; i--) {
			int newEnd=end+i;
			int newStart=newEnd-binSize;
			SingleInterval newPos=new SingleInterval(chr, newStart, newEnd);
			int updatedScore=updateScoreEnd(r, newPos, score);
			int w=binSize;
			double pval=ScanStat.getPValue(updatedScore, lambda, w, genomeLength);
			if(pval<alpha) {
				return newPos;
			}
		}
		return null;*/
	}


	private int getCounts(SingleInterval toAdd, Map<SingleInterval, Integer> endCounts) {
		int rtrn=0;
		for(int i=toAdd.getReferenceStartPosition(); i<=toAdd.getReferenceEndPosition(); i++) {
			SingleInterval pos=new SingleInterval(toAdd.getReferenceName(), i, i);
			rtrn+=get(endCounts,pos);
		}
		return rtrn;
	}

	private int getStartExtension(SingleInterval r, int score) {
		/*String chr=r.getReferenceName();
		int start=r.getReferenceStartPosition();
		SingleInterval newPos=new SingleInterval(chr, start-binSize, start);
		int updatedScore=updateScoreStart(r, newPos, score); 
		
		for(int i=binSize-1; i>0; i--) {
			int newStart=start-i;
			int newEnd=newStart+binSize;
			updatedScore=updatedScore-getCounts(newPos.getReferenceStartPosition(), this.startCounts)+getCounts(newEnd, this.endCounts);	
			newPos=new SingleInterval(chr, newStart, newEnd);
					
			int w=binSize;
			double pval=ScanStat.getPValue(updatedScore, lambda, w, genomeLength);
			if(pval<alpha) {
				return newPos;
			}
		}
		return null;*/
		
		
		int minStart=r.getReferenceStartPosition();
		String chr=r.getReferenceName();
		int start=r.getReferenceStartPosition();
		int updatedScore=score;
		for(int i=1; i<=binSize; i++) {
			int newStart=start-i;
			int newEnd=newStart+binSize;
			updatedScore=updatedScore-getCounts(chr, newEnd+1, this.startCounts)+getCounts(chr, newStart, this.endCounts);
			
			if(updatedScore>cutoff) {minStart=Math.min(minStart, newStart);}
			
			/*int w=binSize;
			double pval=ScanStat.getPValue(updatedScore, lambda, w, genomeLength);
			if(pval<alpha) {
				minStart=Math.min(minStart, newStart);
			}*/
		}
		
		return minStart;
		
		
		/*String chr=r.getReferenceName();
		int start=r.getReferenceStartPosition();
		for(int i=binSize; i>0; i--) {
			int newStart=start-i;
			int newEnd=newStart+binSize;
			SingleInterval newPos=new SingleInterval(chr, newStart, newEnd);
			int updatedScore=updateScoreStart(r, newPos, score);
			int w=binSize;
			double pval=ScanStat.getPValue(updatedScore, lambda, w, genomeLength);
			if(pval<alpha) {
				return newPos;
			}
		}
		return null;*/
	}

	private int getCounts(String chr, int pos, Map<SingleInterval, Integer> startCounts2) {
		SingleInterval newPos=new SingleInterval(chr, pos, pos);
		return getCounts(newPos, startCounts2);
	}


	
	private void write(String save, Collection<SingleInterval> trimmed) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval r: trimmed) {writer.write(r.toShortBED()+"\n");}
		
		writer.close();
	}

	
	private Map<SingleInterval, Integer> significantWindows(File file, int binSize) {
		Map<SingleInterval, Integer> counts=scoreWindows(file, binSize);
		
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		//System.err.println("Total reads "+totalReads);
		this.cutoff=new Double(totalReads).intValue();
		this.lambda=this.totalReads/genomeLength;
		int w=binSize;
		this.numPositionsCovered=positionCount(counts);
		
		Collection<SingleInterval> sigRegions=new TreeSet<SingleInterval>();
		
		for(SingleInterval region: counts.keySet()) {
			double pval=ScanStat.getPValue(counts.get(region), lambda, w, genomeLength);
			if(pval<alpha) {
				this.cutoff=Math.min(cutoff, counts.get(region));
				sigRegions.add(region);
				rtrn.put(region, counts.get(region));
			}
		}
		
		return rtrn;
	}


	


	



	


	private double positionCount(Map<SingleInterval, Integer> counts) {
		double rtrn=0;
		for(SingleInterval r: counts.keySet()) {
			rtrn+=r.size();
		}
		return rtrn;
	}

	
	private Collection<SingleInterval> trim(File file, Collection<SingleInterval> regions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		Map<SingleInterval, int[]> scores=scoreWindows(file, regions);
		
		
		for(SingleInterval region: scores.keySet()) {
			double avg=Statistics.mean(scores.get(region));
			double cutoff=lambda;
			//double cutoff=2*(this.totalReads/this.numPositionsCovered);
			if(avg>0) {
				PoissonDistribution dist=new PoissonDistribution(avg);
				cutoff=Math.max(dist.inverseCumulativeProbability(threshold), cutoff);
			}
			
			
			SingleInterval trimmed=trimEnds(region, scores.get(region), cutoff);
			
			rtrn.add(trimmed);
		}
		
		
		
		return rtrn;
		
	}
	
	
	
	private SingleInterval trim(SingleInterval region, int[] vals, double cutoff) {
		
		double[] norm=normScores(vals, cutoff);
		int[] minmax=MaximumContiguousSubsequence.maxSubSum3(norm);
		
		SingleInterval r=new SingleInterval(region.getReferenceName(), region.getReferenceStartPosition()+minmax[0], region.getReferenceStartPosition()+minmax[1]+1);
		
		return r;
		
		
		/*double[] norm=normScores(vals, cutoff);
		int[] minmax=MaximumContiguousSubsequence.maxSubSum3(norm);
		
		SingleInterval r=new SingleInterval(region.getReferenceName(), region.getReferenceStartPosition()+minmax[0], region.getReferenceStartPosition()+minmax[1]+1);
		
		int len=minmax[1]-minmax[0]+1;
		int[] sub=new int[len];
		
		for(int i=minmax[0]; i<=minmax[1]; i++) {
			sub[i-minmax[0]]=vals[i];
		}
		
		
		double avg=Statistics.mean(sub);
		if(avg>0) {
			PoissonDistribution dist=new PoissonDistribution(avg);
			//System.err.println(region.toUCSC()+" "+cutoff+" "+dist.inverseCumulativeProbability(0.1));
			cutoff=Math.max(dist.inverseCumulativeProbability(threshold), cutoff);
		}
		
		double percent=this.percentLessThan(sub, cutoff);
		
		rtrn.put(r, percent);
		
		return r;*/
	}

	
	
	private SingleInterval trimEnds(SingleInterval region, int[] vals, double cutoff) {
		double[] norm=normScoresBinary(vals, cutoff);
		int startShift=trimStart(norm, binSize);
		int endShift=trimEnd(norm, binSize);
		
		int start=region.getReferenceStartPosition()+startShift;
		int end=region.getReferenceEndPosition()-endShift;
		
		if(end<start) {
			start=region.getReferenceStartPosition();
			end=region.getReferenceEndPosition();
		}
		//System.err.println(region.toUCSC()+" "+vals.length+" "+startShift+" "+endShift);		
		SingleInterval r=new SingleInterval(region.getReferenceName(), start, end);
		
		return r;
		
	}


	private int trimEnd(double[] norm, int binSize2) {
		double[] subset=new double[binSize2];
		
		for(int i=0; i<subset.length; i++) {
			subset[i]=norm[norm.length-i-1];
		}
		
		int endPos=MaximumContiguousSubsequence.maxFromEnd(subset);
		return endPos;
	}

	private int trimStart(double[] norm, int binSize2) {
		double[] subset=new double[binSize2];
		for(int i=0; i<subset.length; i++) {
			subset[i]=norm[i];
		}
		
		return MaximumContiguousSubsequence.maxFromEnd(subset);		
	}

	private Map<SingleInterval, int[]> scoreWindows(File file, Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> trees=makeTree(regions);
		
		
		Map<SingleInterval, int[]> counts=new TreeMap<SingleInterval, int[]>();
		
		for(SingleInterval r: regions) {
			int[] vals=new int[r.size()];
			counts.put(r, vals);
		}
		
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			Collection<SingleInterval> overlappers=getRegions(trees, record);
			for(SingleInterval o: overlappers) {
				increment(counts, o, record);
			}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		reader.close();
		reads.close();
		
		return counts;
		
		
	}
	
	
	


	private void increment(Map<SingleInterval, int[]> counts, SingleInterval region, SAMRecord record) {
		int[] vals=counts.get(region);
		
		int start=Math.max(0, record.getAlignmentStart()-region.getReferenceStartPosition());
		int end=Math.min(region.getGenomicLength(), record.getAlignmentEnd()-region.getReferenceStartPosition());
		
		//System.err.println(region.toUCSC()+" "+record.getAlignmentStart()+" "+record.getAlignmentEnd()+" "+start+" "+end);
		
		for(int i=start; i<end; i++) {
			vals[i]++;
		}
		
	}


	private Map<SingleInterval, Integer> collapse(Map<SingleInterval, Integer> sigRegions) {
		Map<SingleInterval, Integer> rtrn=new TreeMap<SingleInterval, Integer>();
		
		Iterator<SingleInterval> iter=sigRegions.keySet().iterator();
		
		ArrayList<SingleInterval> list=new ArrayList<SingleInterval>();
		int sum=0;
		SingleInterval current=null;
		while(iter.hasNext()) {
			SingleInterval next=iter.next();
			if(ifNext(current, next)) {
				list.add(next);
				sum+=sigRegions.get(next);
				current=next;
			}
			else {
				SingleInterval c=collapseList(list);
				rtrn.put(c, sum);
				list=new ArrayList<SingleInterval>();
				sum=sigRegions.get(next);
				list.add(next);
				current=next;
			}
		}
		
		return rtrn;
	}
	
	
	private Collection<SingleInterval> collapse(Collection<SingleInterval> sigRegions) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Iterator<SingleInterval> iter=sigRegions.iterator();
		
		ArrayList<SingleInterval> list=new ArrayList<SingleInterval>();
		//int sum=0;
		SingleInterval current=null;
		while(iter.hasNext()) {
			SingleInterval next=iter.next();
			if(overlaps(current, next)) {
				list.add(next);
				//sum+=sigRegions.get(next);
				current=next;
			}
			else {
				SingleInterval c=collapseList(list);
				rtrn.add(c);
				list=new ArrayList<SingleInterval>();
				//sum=sigRegions.get(next);
				list.add(next);
				current=next;
			}
		}
		
		return rtrn;
	}
	

	private boolean overlaps(SingleInterval current, SingleInterval next) {
		if(current==null) {return true;}
		
		return current.overlaps(next);
		
		/*if(next.getReferenceStartPosition()==(current.getReferenceEndPosition())) {
			if(next.getReferenceName().equals(current.getReferenceName())) {return true;}
		}
		
		return false;*/
	}
	
	
	private boolean ifNext(SingleInterval current, SingleInterval next) {
		if(current==null) {return true;}
		
		if(next.getReferenceStartPosition()==(current.getReferenceEndPosition())) {
			if(next.getReferenceName().equals(current.getReferenceName())) {return true;}
		}
		
		return false;
	}

	private SingleInterval collapseList(ArrayList<SingleInterval> list) {
		SingleInterval first=list.get(0);
		
		SingleInterval last=list.get(list.size()-1);
		
		return new SingleInterval(first.getReferenceName(), first.getReferenceStartPosition(), last.getReferenceEndPosition());
	}

	
	private double get(Map<SingleInterval, Integer> startCounts2, SingleInterval interval) {
		double count=0;
		if(startCounts2.containsKey(interval)) {count=startCounts2.get(interval);}
		return count;
	}


	
	private Map<String, IntervalTree<SingleInterval>> makeTree(Collection<SingleInterval> regions) {
		Map<String, IntervalTree<SingleInterval>> tree=new TreeMap<String, IntervalTree<SingleInterval>>();
		
		for(SingleInterval r: regions) {
			if(!tree.containsKey(r.getReferenceName())) {
				tree.put(r.getReferenceName(), new IntervalTree<SingleInterval>());
			}
			tree.get(r.getReferenceName()).put(r.getReferenceStartPosition(), r.getReferenceEndPosition(), r);
		}
		
		
		return tree;
	}

	
	

	private double[] normScores(int[] vals, double d) {
		double[] norm=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			norm[i]=(vals[i]-(d));
		}
		
		return norm;
	}
	
	private double[] normScoresBinary(int[] vals, double d) {
		double[] norm=new double[vals.length];
		
		for(int i=0; i<vals.length; i++) {
			double score=(vals[i]-(d));
			if(score>0) {norm[i]=1;}
			else {norm[i]=-1;}
		}
		
		return norm;
	}
	
	private Collection<SingleInterval> getRegions(Map<String, IntervalTree<SingleInterval>> tree, SAMRecord record) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		if(tree.containsKey(record.getReferenceName())) {
			Iterator<SingleInterval> iter=tree.get(record.getReferenceName()).overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(iter.hasNext()) {
				rtrn.add(iter.next());
			}
		}
		
		return rtrn;
	}


	private Map<SingleInterval, Integer> scoreWindows(File file, int binSize) {
		this.startCounts=new TreeMap<SingleInterval, Integer>();
		this.endCounts=new TreeMap<SingleInterval, Integer>();
		
		
		Map<SingleInterval, Integer> counts=new TreeMap<SingleInterval, Integer>();
		SAMFileReader reader=new SAMFileReader(file);
		SAMRecordIterator reads=reader.iterator();
		this.genomeLength=CoordinateSpace.getGenomeLength(reader.getFileHeader());
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			addEnds(record);
			Collection<SingleInterval> allBins=SAMFragment.getSingleInterval(record).allBins(binSize);	
			for(SingleInterval binned: allBins) {
				int score=0;
				if(counts.containsKey(binned)) {score=counts.get(binned);}
				score++;
				counts.put(binned, score);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		this.totalReads=counter;
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


	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			File file=new File(args[0]);
			String save=args[1];
			int binSize= Integer.parseInt(args[2]);
			new SignificantWindows(file, binSize, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=save \n args[2]=binSize";
	
}
