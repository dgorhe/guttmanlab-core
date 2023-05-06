package guttmanlab.core.splicing.speckle.kineticmodel;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class DataByGene {

	public DataByGene(File[] bams, Gene gene, String save) throws IOException {
		List<String> rows=new ArrayList<String>();
		for(int i=0; i<bams.length; i++) {rows.add(bams[i].getName());}
		
		Map<SingleInterval, String> regions=new TreeMap<SingleInterval, String>();
		
		List<String> columns=new ArrayList<String>();
		Iterator<SingleInterval> blocks=gene.getBlocks();
		int i=1;
		while(blocks.hasNext()) {
			SingleInterval exon=blocks.next();
			String name="exon_"+i;
			columns.add(name);
			regions.put(exon, name);
			i++;
		}
		
		Iterator<Annotation> introns=gene.getIntrons().iterator();
		i=1;
		while(introns.hasNext()) {
			SingleInterval intron=introns.next().getSingleInterval();
			String name="intron_"+i;
			columns.add(name);
			columns.add("spliced_"+i);
			regions.put(intron, name);
			i++;
		}
		
		
		MatrixWithHeaders matrix=new MatrixWithHeaders(rows, columns);
		
		
		for(int index=0; index<bams.length; index++) {
			System.err.println(bams[index].getName()+" "+(index+1)+" "+bams.length);
			Map<String, Double> counts=count(bams[index], gene, regions);
			
			for(String col: counts.keySet()) {
				matrix.set(bams[index].getName(), col, counts.get(col));
			}
			
		}
		
		matrix.write(save);
	}
	
	
	
	private Map<String, Double> count(File bam, Gene gene, Map<SingleInterval, String> regions) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		List<SAMRecord> list=new ArrayList<SAMRecord>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			if(!record.getDuplicateReadFlag()) {
			
				
				
				if(record.getReferenceName().equals(gene.getReferenceName())) {
					if(record.getAlignmentStart()>=gene.getReferenceStartPosition() && record.getAlignmentEnd()<=gene.getReferenceEndPosition()) {
						if(SAMFragment.getOrientation(record).equals(gene.getOrientation())) {
							list.add(record);
						}
					}
				}
				
				
				counter++;
			}
				
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		double total=counter;
		reader.close();
		reads.close();
		
		for(SAMRecord read: list) {
			SAMFragment frag=new SAMFragment(read);
			if(frag.isSpliced()) {
				Collection<SingleInterval> readJunctions=getJunctions(frag);
				for(SingleInterval region: regions.keySet()) {
					String name1=regions.get(region);
					if(name1.startsWith("intron")) {
						String name="spliced_"+name1.split("_")[1];
						if(readJunctions.contains(region)) {
							increment(name, rtrn);
						}
					}
				}
			}
			else {
				for(SingleInterval region: regions.keySet()) {
					String name=regions.get(region);
					if(overlap(read, region)) {
						increment(name, rtrn);
					}
				}
			}
		}
		
		
		Map<String, Double> norm=new TreeMap<String, Double>();
		for(String name: rtrn.keySet()) {
			double count=rtrn.get(name);
			double val=(count/total)*1000000;
			norm.put(name, val);
		}
		
		
		
	
		return norm;
	}

	
	


	private boolean overlap(SAMRecord read, SingleInterval region) {
		if(read.getAlignmentStart()>=region.getReferenceStartPosition() && read.getAlignmentEnd()<=region.getReferenceEndPosition()) {return true;}
		return false;
	}



	private void increment(String name, Map<String, Integer> rtrn) {
		int counter=0;
		if(rtrn.containsKey(name)) {counter=rtrn.get(name);}
		counter++;
		rtrn.put(name, counter);
		
	}



	private static int[] assign(Gene g, SAMFragment read) {
		
		int[] rtrn=new int[3];
		
		String state=assignState(read, g);
		if(state.equals("splicedJunction")) {
			rtrn[0]++;
		}
		if(state.equals("splicedNoJunction")) {rtrn[1]++;}
		if(state.equals("unspliced")) {
			rtrn[2]++;
		}
		//Pair<Integer> rtrn=new Pair<Integer>(exonCount, intronCount);
		
		return rtrn;
		
		
	}
	
	
	private static String assignState(SAMFragment read, Gene gene) {
		//If read is spliced, check if splice junction matches gene junction, if so --> spliced
		if(read.isSpliced()) {
			if(junctionMatch(read, gene)) {return "splicedJunction";}
			else {return "amb";}
		}
		
		for(Annotation intron: gene.getIntrons()) {
			if(read.overlaps(intron)) {
				//TODO Trim read by 3 bases
				return "unspliced";
			}
		}
		//return "spliced";
		return "splicedNoJunction";
	}

	private static boolean junctionMatch(SAMFragment read, Gene gene) {
		if(read.getOrientation().equals(gene.getOrientation())) {
			Collection<SingleInterval> readJunctions=getJunctions(read);
			Collection<SingleInterval> geneJunctions=getJunctions(gene);
			for(SingleInterval readJunction: readJunctions) {
				if(geneJunctions.contains(readJunction)){return true;}
			}	
		}
		return false;
	}
	
	private static Collection<SingleInterval> getJunctions(Annotation gene) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		Collection<Annotation> blocks=gene.getIntrons();
		
		for(Annotation intron: blocks) {
			SingleInterval intronSI=new SingleInterval(intron.getReferenceName(), intron.getReferenceStartPosition(), intron.getReferenceEndPosition());
			intronSI.setOrientation(gene.getOrientation());
			rtrn.add(intronSI);
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
		
			Gene transcript=BEDFileIO.loadRegionsFromFile(args[0]).iterator().next();
			File[] files=new File(args[1]).listFiles();
			new DataByGene(files, transcript, args[2]);
		
		}else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=bed \n args[1]=bam files \n args[2]=save";
}
