package guttmanlab.core.sars;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SplicingAnalysisPerGene {
	
	//TODO Add excluded regions

	public SplicingAnalysisPerGene(File bam, Collection<Gene> genes, String save) throws IOException{
		Map<String, IntervalTree<Annotation>> exons=getExonTree(genes);
		Map<String, IntervalTree<Annotation>> introns=getIntronTree(genes);
		
		Map<Annotation, Integer> intronCounts=new TreeMap<Annotation, Integer>();
		Map<Annotation, Integer> exonCounts=new TreeMap<Annotation, Integer>();
		
		SAMFileReader inputReader= new SAMFileReader(bam);
		/*SAMFileWriter ambWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(inputReader.getFileHeader(), false, new File(save+".amb.bam"));
		SAMFileWriter exonWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(inputReader.getFileHeader(), false, new File(save+".exon.bam"));
		SAMFileWriter intronWriter=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(inputReader.getFileHeader(), false, new File(save+".intron.bam"));
		*/
		
		double intronCount=0;
		double exonCount=0;
		
		int totalCount=0;
		SAMRecordIterator reads=inputReader.iterator();
		while(reads.hasNext()){
			SAMRecord read=reads.next();
			boolean exon=overlaps(exons, read);
			boolean intron=overlaps(introns, read);
			
			if(exon && intron){
				
			}
			else if(exon){
				Collection<Annotation> regions=overlappingList(exons, read);
				for(Annotation region: regions){
					int count=0;
					if(exonCounts.containsKey(region)){
						count=exonCounts.get(region);
					}
					count++;
					exonCounts.put(region, count);
				}
				exonCount++;
			}
			else if(intron){
				Collection<Annotation> regions=overlappingList(introns, read);
				//System.err.println(regions.size());
				for(Annotation region: regions){
					int count=0;
					if(intronCounts.containsKey(region)){
						count=intronCounts.get(region);
					}
					count++;
					intronCounts.put(region, count);
				}
				
				intronCount++;
			}
			else{System.err.println("Amb");}
			totalCount++;
			if(totalCount%100000 ==0){System.err.println(totalCount);}
		}
				
				
		reads.close();
		inputReader.close();
		
		write(save+".intron.counts", intronCounts);
		write(save+".exon.counts", exonCounts);
		
		double exonRatio=10000*(exonCount/(double)totalCount);
		double IntronRatio=10000*(intronCount/(double)totalCount);
		
		//System.err.println(exonRatio+" "+IntronRatio);
	}

	private Collection<Annotation> overlappingList(Map<String, IntervalTree<Annotation>> introns, SAMRecord read) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		String chr=read.getReferenceName();
		if(!chr.startsWith("chr")){chr="chr"+read.getReferenceName();}
		System.err.println(chr+" "+read.getReferenceName());
		if(introns.containsKey(chr)){
			IntervalTree<Annotation> tree=introns.get(chr);
			Iterator<Annotation> iter=tree.overlappingValueIterator(read.getAlignmentStart(), read.getAlignmentEnd());
			while(iter.hasNext()){rtrn.add(iter.next());}
			System.err.println(chr+" "+read.getAlignmentStart()+" "+ read.getAlignmentEnd()+" "+rtrn.size());
		}
		
		
		return rtrn;
	}

	private void write(String save, Map<Annotation, Integer> intronCounts) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation region: intronCounts.keySet()){
			writer.write(region.tobedgraph(intronCounts.get(region))+"\t"+region.size()+"\n");
		}
		
		writer.close();
	}

	private Map<String, IntervalTree<Annotation>> getExonTree(Collection<Gene> genes) {
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		for(Gene gene: genes){
			//System.err.println("exons "+gene.getName());
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			if(rtrn.containsKey(gene.getReferenceName())){tree=rtrn.get(gene.getReferenceName());}
			Iterator<SingleInterval> exons=gene.getBlocks();
			while(exons.hasNext()){
				SingleInterval exon=exons.next();
				tree.put(exon.getReferenceStartPosition(),  exon.getReferenceEndPosition(), exon);
			}
			rtrn.put(gene.getReferenceName(), tree);
		}
		
		return rtrn;
	}
	
	private Map<String, IntervalTree<Annotation>> getIntronTree(Collection<Gene> genes) {
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		for(Gene gene: genes){
			//System.err.println(gene.getName());
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			if(rtrn.containsKey(gene.getReferenceName())){tree=rtrn.get(gene.getReferenceName());}
			Collection<Annotation> introns=gene.getIntrons();
			for(Annotation exon:introns){
				tree.put(exon.getReferenceStartPosition(),  exon.getReferenceEndPosition(), exon);
			}
			rtrn.put(gene.getReferenceName(), tree);
		}
		
		return rtrn;
	}

	private boolean overlaps(Map<String, IntervalTree<Annotation>> exons, SAMRecord read) {
		String chr="chr"+read.getReferenceName();
		if(exons.containsKey(chr)){
			IntervalTree<Annotation> tree=exons.get(chr);
			return tree.hasOverlappers(read.getAlignmentStart(), read.getAlignmentEnd());
		}
		return false;
	}

	private void write(SAMRecord read, SAMFileWriter ambWriter) {
		ambWriter.addAlignment(read);
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File bam=new File(args[0]);
			Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[1]);
			String save=args[2];
			new SplicingAnalysisPerGene(bam, genes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam \n args[1]=gene bed file \n args[2]=save";
}
