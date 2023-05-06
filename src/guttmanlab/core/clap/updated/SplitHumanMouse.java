package guttmanlab.core.clap.updated;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplitHumanMouse {

	private static void splitHumanMouse(File bam, String save, Map<String, IntervalTree<Gene>> trees) {
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		SAMFileWriter writerHuman=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save+".human.bam"));
		SAMFileWriter writerMouse=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(reader.getFileHeader(), false, new File(save+".mouse.bam"));
		
		int humanCount=0;
		int mouseCount=0;
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(!record.getDuplicateReadFlag()) {
				String chr=record.getReferenceName();
				String trimmed=trim(chr);
				//record.setReferenceName(trimmed);
				if(chr.contains("_human")) {
					if(CoordinateSpace.HG38.getRefSizes().containsKey(trimmed)) {
						boolean hasGene=overlappingGene(record, trees);
						if(hasGene) {
							writerHuman.addAlignment(record);	
							humanCount++;
						}
					}
				}
				if(chr.contains("_mouse")) {
					if(CoordinateSpace.MM10.getRefSizes().containsKey(trimmed)) {
						boolean hasGene=overlappingGene(record, trees);
						if(hasGene) {
							writerMouse.addAlignment(record);
							mouseCount++;
						}
					}
				}
				
				counter++;
				if(counter%1000000==0) {System.err.println(counter);}
			}
		}
		
		
		
		reader.close();
		reads.close();
		writerHuman.close();
		writerMouse.close();
		
		//int totalBinNumber=score(bam, 100).size();
		int humanBinNumber=score( new File(save+".human.bam"), 100).size();
		int mouseBinNumber=score( new File(save+".mouse.bam"), 100).size();
		
		System.out.println(bam.getName()+"\t"+counter+"\t"+humanCount+"\t"+mouseCount+"\t"+humanBinNumber+"\t"+mouseBinNumber);	
	}
	
	private static boolean overlappingGene(SAMRecord record, Map<String, IntervalTree<Gene>> trees) {
		SAMFragment fragment=new SAMFragment(record);
		String chr=record.getReferenceName();
		if(trees.containsKey(chr)) {
			IntervalTree<Gene> tree=trees.get(chr);
			Iterator<Gene> iter=tree.overlappingValueIterator(record.getAlignmentStart(), record.getAlignmentEnd());
			while(iter.hasNext()){
				Gene g=iter.next();
				SingleInterval gSI=g.getGenomicRegion();
				if(gSI.overlaps(fragment) || fragment.overlaps(gSI)) {return true;}
			}
		}
		return false;
	}

	private static void writeBins(File bam, int binSize) {
		Map<SingleInterval, Integer> map=score(bam, binSize);
		
		for(SingleInterval region: map.keySet()) {
			System.out.println(region.toBedgraph(map.get(region)));
		}
		
	}
	
	
	private static void scoreChromosomes(File bam, String save) throws IOException {
		Map<String, Integer> map=scoreChr(bam);
		
		FileWriter writer=new FileWriter(save);
		
		for(String chr: map.keySet()) {
			writer.write(chr+"\t"+map.get(chr)+"\n");
		}
		
		writer.close();
	}
	

	private static Map<SingleInterval, Integer> score(File bam, int binSize){
		SAMFileReader inputReader= new SAMFileReader(bam);
		Map<SingleInterval, Integer> positionCount=new TreeMap<SingleInterval, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getReadUnmappedFlag() && !read.getDuplicateReadFlag()){
				//if(!read.getProperPairFlag()) {System.err.println(read.getReadName());}
				Collection<SingleInterval> allBins=SAMFragment.allBins(read, binSize);
				for(SingleInterval bin: allBins) {
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
		
		//System.out.println(bam.getName()+"\t"+totalCount+"\t"+positionCount.keySet());
		
		return positionCount;
	}
	
	
	private static Map<String, Integer> scoreChr(File bam){
		SAMFileReader inputReader= new SAMFileReader(bam);
		Map<String, Integer> positionCount=new TreeMap<String, Integer>();
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(!read.getDuplicateReadFlag()){
				String chr=read.getReferenceName();
				int score=0;
				if(positionCount.containsKey(chr)){score=positionCount.get(chr);}
				score++;
				positionCount.put(chr, score);
			}
			
			totalCount++;
			if(totalCount%1000000 ==0){System.err.println(totalCount+" "+read.getReferenceName());}
		}
				
		reads.close();
		inputReader.close();
		
		return positionCount;
	}
	
	private static void filterBins(String string, int parseInt) throws NumberFormatException, IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(string)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			int score=Integer.parseInt(tokens[3]);
			if(score>parseInt) {System.out.println(nextLine);}
		}
		reader.close();
		
	}
	
	private static List<File> getBams(File[] listFiles) {
		List<File> rtrn=new ArrayList<File>();
		for(int i=0; i<listFiles.length; i++) {
			File file=listFiles[i];
			if(file.getName().endsWith(".bam")) {
				rtrn.add(file);
			}
		}
		return rtrn;
	}

	
	
	private static String trim(String chr) {
		String rtrn=chr.replace("_human", "");
		rtrn=rtrn.replace("_mouse", "");
		return rtrn;
	}
	
	
	private static void scoreChromosomes(List<File> bams, String save) throws IOException {
		Map<String, Integer>[] counts=new Map[bams.size()];
		
		int i=0;
		for(File bam: bams) {
			System.err.println(bam.getAbsolutePath());
			counts[i]=scoreChr(bam);
			i++;
		}
		
		write(save, counts, bams);
		
	}
	
	private static void write(String save, Map<String, Integer>[] counts, List<File> bams) throws IOException {
		Collection<String> allChr=new TreeSet<String>();
		for(int i=0; i<counts.length; i++) {allChr.addAll(counts[i].keySet());}
		
		FileWriter writer=new FileWriter(save);
		writer.write("chromosome");
		for(File bam: bams) {writer.write("\t"+bam.getName());}
		writer.write("\n");
		
		
		int[] human=sumHuman(counts, allChr);
		int[] mouse=sumMouse(counts, allChr);
		int[] repeat=sumRepeat(counts, allChr);
		int[] unmapped=sumUnmapped(counts, allChr);
		
		writer.write("human");
		for(int i=0; i<human.length; i++) {
			writer.write("\t"+human[i]);
		}
		writer.write("\n");
		
		writer.write("mouse");
		for(int i=0; i<mouse.length; i++) {
			writer.write("\t"+mouse[i]);
		}
		writer.write("\n");
		
		writer.write("repeat");
		for(int i=0; i<repeat.length; i++) {
			writer.write("\t"+repeat[i]);
		}
		writer.write("\n");
		
		writer.write("unmapped");
		for(int i=0; i<unmapped.length; i++) {
			writer.write("\t"+unmapped[i]);
		}
		writer.write("\n");
		
		/*for(String chr: allChr) {
			writer.write(chr);
			for(int i=0; i<counts.length; i++) {
				int count=getCount(counts[i], chr);
				writer.write("\t"+count);
			}
			writer.write("\n");
		}*/
		
		writer.close();
	}

	private static int[] sumRepeat(Map<String, Integer>[] counts, Collection<String> allChr) {
		int[] rtrn=new int[counts.length];
		
		for(String chr: allChr) {
			if(!chr.endsWith("_human") && !chr.endsWith("*") && !chr.endsWith("_mouse")) {
				for(int i=0; i<counts.length; i++) {
					rtrn[i]+=getCount(counts[i], chr);
				}
			}
		}
		return rtrn;
	}

	private static int[] sumUnmapped(Map<String, Integer>[] counts, Collection<String> allChr) {
		int[] rtrn=new int[counts.length];
		
		for(String chr: allChr) {
			if(chr.endsWith("*")) {
				for(int i=0; i<counts.length; i++) {
					rtrn[i]+=getCount(counts[i], chr);
				}
			}
		}
		return rtrn;
	}

	private static int[] sumHuman(Map<String, Integer>[] counts, Collection<String> allChr) {
		int[] rtrn=new int[counts.length];
		
		for(String chr: allChr) {
			if(chr.endsWith("_human")) {
				for(int i=0; i<counts.length; i++) {
					rtrn[i]+=getCount(counts[i], chr);
				}
			}
		}
		return rtrn;
	}
	
	private static int[] sumMouse(Map<String, Integer>[] counts, Collection<String> allChr) {
		int[] rtrn=new int[counts.length];
		
		for(String chr: allChr) {
			if(chr.endsWith("_mouse")) {
				for(int i=0; i<counts.length; i++) {
					rtrn[i]+=getCount(counts[i], chr);
				}
			}
		}
		return rtrn;
	}

	private static int getCount(Map<String, Integer> map, String chr) {
		if(map.containsKey(chr)) {return map.get(chr);}
		return 0;
	}

	public static void main(String[] args) throws NumberFormatException, IOException {
		if(args.length>1) {
		File bam=new File(args[0]);
		String save=args[1];
		scoreChromosomes(bam, save);
		}
		else {System.err.println(usage);}
	}

	static String usage=" args[0]=bam \n args[1]=save";

	
	
	
}
