package guttmanlab.core.trip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.Sequence;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class PullBarcodeSequence {
	
	static int maxDistance=2;
	static int maxMismatch=5;
	private static int collapseDistance=1000000;

	/*private static Map<String, Integer> pullBarcodes(File fastq2) {
		Map<String, Integer> barcodeCount=new TreeMap<String, Integer>();
		FastqReader f=new FastqReader(fastq2);
			
		Iterator<FastqRecord> iter=f.iterator();
			
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String seqBases=seq.getReadString();
			String barcode=seqBases.substring(0, 16);
			int count=0;
			if(barcodeCount.containsKey(barcode)) {count=barcodeCount.get(barcode);}
			count++;
			barcodeCount.put(barcode, count);
			counter++;
			if(counter%1000000==0) {System.err.println(counter+"\t"+barcode+"\t"+barcode.length());}
		}
			
		f.close();
		return barcodeCount;
		//write(barcodeCount);
		
	}*/
	
	private static Map<String, Collection<String>> pullBarcodes(File fastq2) {
		String commonSeq="CTGCAGTTA";
		
		Map<String, Collection<String>> barcodeCount=new TreeMap<String, Collection<String>>();
		FastqReader f=new FastqReader(fastq2);
			
		Iterator<FastqRecord> iter=f.iterator();
			
		int counter=0;
		while(iter.hasNext()) {
			FastqRecord seq=iter.next();
			String seqBases=seq.getReadString();
			String barcode=seqBases;
			if(!barcodeCount.containsKey(barcode)) {barcodeCount.put(barcode, new TreeSet<String>());}
			Collection<String> list=barcodeCount.get(barcode);
			String name=seq.getReadHeader().split(" ")[0];
			list.add(name);
			barcodeCount.put(barcode, list);
			
			/*int commonPos=seqBases.indexOf(commonSeq);
			
			if(commonPos==37) {
				String barcode=seqBases.substring(commonPos-16, commonPos);
				//System.err.println(commonPos+" "+barcode+" "+barcode.length());	
				if(!barcodeCount.containsKey(barcode)) {barcodeCount.put(barcode, new TreeSet<String>());}
				Collection<String> list=barcodeCount.get(barcode);
				list.add(seq.getReadHeader());
				barcodeCount.put(barcode, list);
			}
			else {}*/
			
			counter++;
			if(counter%1000000==0) {System.err.println(counter);}
		}
		
		f.close();
		return barcodeCount;
			
		
		
	}
	
	private static void write(Map<String, Integer> barcodeCount) {
		for(String barcode: barcodeCount.keySet()) {
			int count=barcodeCount.get(barcode);
			System.out.println(barcode+"\t"+count);
		}
		
	}
	
	private static void writeList(Map<String, Collection<String>> cDNABarcode) {
		for(String barcode: cDNABarcode.keySet()) {
			System.out.println(barcode+"\t"+cDNABarcode.get(barcode).size());
		}
		
	}
	
	/*private static void merge(Map<String, Integer> map1, Map<String, Integer> map2) {
		Collection<String> allBarcodes=new TreeSet<String>();
		allBarcodes.addAll(map1.keySet());
		allBarcodes.addAll(map2.keySet());
		
		int count=0;
		for(String barcode: allBarcodes) {
			int val1=get(map1, barcode);
			int val2=get(map2, barcode);
			if(val2==0) {val2=getClosest(map1, barcode);}
			if(val1!=0 && val2!=0) {count++;}
			System.out.println(barcode+"\t"+val1+"\t"+val2);
		}
		
		System.err.println(count);
	}*/

	private static int getClosest(Map<String, Integer> map1, String barcode) {
		int minDistance=2;
		for(String barcode2: map1.keySet()) {
			if(distance(barcode2, barcode)<minDistance) {
				//System.err.println(barcode+" "+barcode2);
				return map1.get(barcode2);
			}
		}
		return 0;
	}
	
	private static String getClosest(Set<String> map1, String barcode) {
		int minDistance=2;
		for(String barcode2: map1) {
			if(distance(barcode2, barcode)<minDistance) {
				//System.err.println(barcode+" "+barcode2);
				return barcode2;
			}
		}
		return null;
	}

	private static int distance(String seq1, String seq2) {
		char[] char1=seq1.toCharArray();
		char[] char2=seq2.toCharArray();
		
		int distance=0;
		for(int i=0; i<char1.length; i++) {
			if(char1[i]!=char2[i]) {distance++;}
		}
		return distance;
	}

	private static int get(Map<String, Integer> map1, String barcode) {
		int rtrn=0;
		
		if(map1.containsKey(barcode)) {rtrn=map1.get(barcode);}
		
		return rtrn;
	}
	
	private static Map<String, Pair<Integer>> getSpliceState(File splicedBam, File unsplicedBam) {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		int counter=0;
		
		SAMFileReader inputReader= new SAMFileReader(splicedBam);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			if(read.getAlignmentStart()<93 && read.getAlignmentEnd()>93) {
				String name=read.getReadName();
				Pair<Integer> pair=new Pair<Integer>(1,0);
				rtrn.put(name, pair);
			}
			else {System.err.println(SAMFragment.getSingleInterval(read).toUCSC());}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
				
		reads.close();
		inputReader.close();
		
		
		counter=0;
		inputReader= new SAMFileReader(unsplicedBam);
		reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			Pair<Integer> pair=new Pair<Integer>(0,1);
			rtrn.put(name, pair);
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
				
		reads.close();
		inputReader.close();
		
		
		
		
		return rtrn;
	}
	
	private static Map<String, SingleInterval> getReadToLocation(File bam) {
		Map<String, SingleInterval> rtrn=new TreeMap<String, SingleInterval>();
		int counter=0;
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			String name=read.getReadName();
			if(Integer.parseInt(read.getAttribute("nM").toString())<maxMismatch) {
			
				SingleInterval r=SAMFragment.getSingleInterval(read).bin(collapseDistance);
				rtrn.put(name, r);
			}
			
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
				
		reads.close();
		inputReader.close();
		
		
		
		
		
		
		
		return rtrn;
	}
	
	
	private static Map<String, Pair<Integer>> merge(Map<String, Collection<String>> cDNABarcode, Map<String, Pair<Integer>> spliceState) {
		Map<String, Pair<Integer>> rtrn=new TreeMap<String, Pair<Integer>>();
		for(String barcode: cDNABarcode.keySet()) {
			Collection<String> readNames=cDNABarcode.get(barcode);
			Pair<Integer> sum=sum(spliceState, readNames);
			//System.out.println(barcode+"\t"+sum.getValue1()+"\t"+sum.getValue2());
			rtrn.put(barcode, sum);
		}
		return rtrn;
	}
	
	private static Map<String, Collection<SingleInterval>> mergeGenomic(Map<String, Collection<String>> cDNABarcode, Map<String, SingleInterval> genomicRegions) {
		int countUnmapped=0;
		int countUnique=0;
		int countLarge=0;
		
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		for(String barcode: cDNABarcode.keySet()) {
			Collection<String> readNames=cDNABarcode.get(barcode);
			Collection<SingleInterval> regions=getRegions(genomicRegions, readNames);
			rtrn.put(barcode, regions);
			//Pair<Integer> sum=sum(spliceState, readNames);
			if(regions.size()==1) {countUnique++;}
			if(regions.size()>1) {countLarge++;}
			if(regions.isEmpty()) {countUnmapped++;}
			
			if(!regions.isEmpty()) {
				System.out.print(barcode+"\t"+regions.size());
				for(SingleInterval r: regions) {System.out.print("\t"+r.toUCSC());}
				System.out.println();
			}
		}
		
		System.err.println(countUnmapped+" "+countUnique+" "+countLarge);
		
		return rtrn;
	}

	

	

	private static Collection<SingleInterval> getRegions(Map<String, SingleInterval> genomicRegions, Collection<String> readNames) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(String read: readNames) {
			if(genomicRegions.containsKey(read)) {
				rtrn.add(genomicRegions.get(read));
			}
		}
		
		return rtrn;
	}

	private static Pair<Integer> sum(Map<String, Pair<Integer>> spliceState, Collection<String> readNames) {
		Pair<Integer> rtrn=new Pair<Integer>(0,0);
		for(String name: readNames) {
			//System.err.println(name);
			if(spliceState.containsKey(name)) {
				Pair<Integer> val=spliceState.get(name);
				rtrn.setValue1(rtrn.getValue1()+val.getValue1());
				rtrn.setValue2(rtrn.getValue2()+val.getValue2());
			}
		}
		return rtrn;
	}
	
	private static void mergeGenomicAndSpliced(Map<String, Collection<SingleInterval>> barcodeToGenomic, Map<String, Pair<Integer>> barcodeToSpliceCounts, String save) throws IOException {
		barcodeToGenomic=filter(barcodeToGenomic);
		
		
		Map<String, String> barcodeMatches=match(barcodeToSpliceCounts.keySet(), barcodeToGenomic.keySet());
		FileWriter writer=new FileWriter(save);
		int counter=0;
		for(String barcode: barcodeMatches.keySet()) {
			String barcode2=barcodeMatches.get(barcode);
			if(barcode2!=null) {
				Collection<SingleInterval> regions=barcodeToGenomic.get(barcode2);
				int distance=distance(barcode, barcode2);
				Pair<Integer> splicedNums=barcodeToSpliceCounts.get(barcode);
				write(writer, barcode, barcode2, distance, regions, splicedNums);
			}
			counter++;
		}
		if(counter%1000==0) {System.err.println(counter);}
		writer.close();
	}

	private static Map<String, Collection<SingleInterval>> filter(Map<String, Collection<SingleInterval>> barcodeToGenomic) {
		Map<String, Collection<SingleInterval>> rtrn=new TreeMap<String, Collection<SingleInterval>>();
		
		for(String b: barcodeToGenomic.keySet()) {
			Collection<SingleInterval> regions=barcodeToGenomic.get(b);
			if(regions.size()==1) {rtrn.put(b, barcodeToGenomic.get(b));}
		}
		
		System.err.println("filtered "+barcodeToGenomic.size()+" "+rtrn.size());
		
		return rtrn;
	}

	private static String getClosestBarcode(String barcode, Map<String, Collection<SingleInterval>> barcodeToGenomic) {
		if(barcodeToGenomic.containsKey(barcode)) {return barcode;}
		
		String barcode2=get(barcodeToGenomic, barcode, maxDistance);
		if(barcode2!=null) {return barcode2;}
		
		
		return null;
	}

	

	private static String get(Map<String, Collection<SingleInterval>> barcodeToGenomic, String barcode, int distance) {
		for(String barcode2: barcodeToGenomic.keySet()) {
			if(distance(barcode2, barcode)<=distance) {
				return barcode2;
			}
		}
		return null;
	}

	private static void write(FileWriter writer, String barcode, String barcode2, int distance,	Collection<SingleInterval> regions, Pair<Integer> splicedNums) throws IOException {
		writer.write(barcode+"\t"+barcode2+"\t"+distance+"\t"+splicedNums.getValue1()+"\t"+splicedNums.getValue2());
		for(SingleInterval r: regions) {
			writer.write("\t"+r.toUCSC());
		}
		writer.write("\n");
	}

	
	private static Map<String, String> match(Set<String> barcodes1, Set<String> barcodes2) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		
		for(String barcode1: barcodes1) {
			if(barcodes2.contains(barcode1)) {
				rtrn.put(barcode1, barcode1);
				//System.out.println(barcode1+"\t"+barcode1+"\t"+0);
			}
			else {
				Collection<String> mismatches=mismatch(barcode1);
				String barcode2=get(mismatches, barcodes2);
				if(barcode2!=null) {rtrn.put(barcode1, barcode2);}
				else {
					for(String m: mismatches) {
						Collection<String> mismatches2= mismatch(m);
						barcode2=get(mismatches2, barcodes2);
						if(barcode2!=null) {rtrn.put(barcode1, barcode2);}
					}
				}
				/*String barcode2=getClosest(barcodes2, barcode1);
				if(barcode2!=null) {
					int distance=distance(barcode1, barcode2);
					System.out.println(barcode1+"\t"+barcode2+"\t"+distance);
				}*/
			}
		}
		return rtrn;
		//System.err.println(unmatched.size());
	}
	
	
	private static String get(Collection<String> mismatches, Set<String> barcodes2) {
		for(String m: mismatches) {
			if(barcodes2.contains(m)) {
				return m;
			}
		}
		return null;
	}

	private static Collection<String> mismatch(String seq) {
		Collection<String> rtrn=new TreeSet<String>();
		for(int i=0; i<seq.length(); i++){
			String first=seq.substring(0, i);
			String second=seq.substring(i+1, seq.length());
			rtrn.add(first+"A"+second);
			rtrn.add(first+"C"+second);
			rtrn.add(first+"G"+second);
			rtrn.add(first+"T"+second);
			rtrn.add(first+"N"+second);
		}
		return rtrn;
	}
	
	
	private static Map<String, Collection<String>> closeBarcodes(Set<String> barcodes) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String barcode1: barcodes) {
			Collection<String> mismatches=mismatch(barcode1);
			Collection<String> list=new TreeSet<String>();
			list.add(barcode1);
			for(String barcode2: mismatches) {
				if(barcodes.contains(barcode2)) {list.add(barcode2);}
			}
			rtrn.put(barcode1, list);
			//System.out.println(barcode1+"\t"+list.size());
		}
		
		
		List<Collection<String>> sets=collapse(rtrn);
		
		for(Collection<String> set: sets) {
			String line=""+set.size();
			for(String b: set) {line+="\t"+b;}
			System.out.println(line);
		}
		
		/*Collection<String> uniqueSet=new TreeSet<String>();
		for(String k: rtrn.keySet()) {
			String unique=toString(rtrn.get(k));
			uniqueSet.add(unique);
		}
		
		System.err.println(uniqueSet.size());
		
		for(String u: uniqueSet) {
			System.out.println(u);
		}*/
		
		return rtrn;
	}

	
	
	private static List<Collection<String>> collapse(Map<String, Collection<String>> map) {
		//go through each kmer and get close hits and close hits of those and collapse
		
		List<Collection<String>> sets=new ArrayList<Collection<String>>();
		for(String barcode: map.keySet()) {
			Collection<String> connected=map.get(barcode);
			sets.add(connected);
		}
		return sets;
	}

	



	private static Map<String, Pair<Integer>> collapseByGenomicPosition(String file) throws IOException {
		Map<String, Pair<Integer>> vals=new TreeMap<String, Pair<Integer>>();
		List<String> lines=BEDFileIO.loadLines(file);
		for(String line: lines) {
			String[] tokens=line.split("\t");
			String position=tokens[5];
			Pair<Integer> val=new Pair<Integer>();
			val.setValue1(Integer.parseInt(tokens[3]));
			val.setValue2(Integer.parseInt(tokens[4]));
			if(vals.containsKey(position)) {
				Pair<Integer> other=vals.get(position);
				val.setValue1(val.getValue1()+other.getValue1());
				val.setValue2(val.getValue2()+other.getValue2());
			}
			vals.put(position, val);
		}
		return vals;
	}


	private static void print(Map<String, Pair<Integer>> barcodeToSpliceCounts) {
		for(String barcode: barcodeToSpliceCounts.keySet()) {
			Pair<Integer> pair=barcodeToSpliceCounts.get(barcode);
			System.out.println(barcode+"\t"+pair.getValue1()+"\t"+pair.getValue2());
		}
		
	}

	
	public static void main(String[] args) throws IOException {
		//Map<String, Collection<String>> cDNABarcode=pullBarcodes(new File(args[0])); //barcode to names
		
		//writeList(cDNABarcode);
		
		//Map<String, Collection<String>> gDNABarcode=pullBarcodes(genomicFile); //barcode to names
		
		
		
		if(args.length>2) {
			//Splicing
			File cDNAFile=new File(args[0]);
			Map<String, Collection<String>> cDNABarcode=pullBarcodes(cDNAFile); //barcode to names
			Map<String, Pair<Integer>> spliceState=getSpliceState(new File(args[1]), new File(args[2])); //name to splice state
			Map<String, Pair<Integer>> barcodeToSpliceCounts=merge(cDNABarcode, spliceState);
			
			
			print(barcodeToSpliceCounts);
			
			//genomic DNA
			/*File r2=new File(args[3]);
			Map<String, Collection<String>> gDNABarcode=pullBarcodes(r2); //barcode to names
			Map<String, SingleInterval> genomicLocation=getReadToLocation(new File(args[4]));
			Map<String, Collection<SingleInterval>> barcodeToGenomic=mergeGenomic(gDNABarcode, genomicLocation);
		
			String save=args[5];
			
			mergeGenomicAndSpliced(barcodeToGenomic, barcodeToSpliceCounts, save);*/
			
			
			
			
		
			
			//TODO match within 2 edit distance
			//TODO collapse cDNAs within edit distance
			//TODO filter unmapped genomic barcodes
			//TODO compute hamming distance of cDNA and gDNA barcodes
			
			
			
			
		}
		else {System.err.println(usage);}
		
	}



	



	static String usage=" args[0]=cDNA R2 fastq \n args[1]=spliced alignments (bam) \n args[2]=unspliced alignments (bam) \n args[3]=gDNA R2 (fastq) \n args[4]=gDNA alignments (bam) \n args[5]=save";

	

	
	
}
