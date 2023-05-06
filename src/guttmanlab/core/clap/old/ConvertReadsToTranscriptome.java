package guttmanlab.core.clap.old;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.sequence.Sequence;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;

public class ConvertReadsToTranscriptome {

	public ConvertReadsToTranscriptome(File bamFile, AnnotationCollection<Gene> genes, File save){
		BAMPairedFragmentCollection bam=new BAMPairedFragmentCollection(bamFile);
		CloseableIterator<PairedMappedFragment<SAMFragment>> reads=bam.sortedIterator();
		//SAMFileHeader header=new SAMFileHeader();
		SAMFileHeader header=createHeader(genes);
		SAMFileWriter writer=new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, save);
		
		while(reads.hasNext()){
			PairedMappedFragment<SAMFragment> read=reads.next();
			Collection<Pair<SAMRecord>> newReads=convert(read, genes, header);
			for(Pair<SAMRecord> newRead: newReads){
				SAMRecord read1=newRead.getValue1();
				SAMRecord read2=newRead.getValue2();
				writer.addAlignment(read1);
				writer.addAlignment(read2);
			}
		}
		
		reads.close();
		writer.close();
	}

	private SAMFileHeader createHeader(AnnotationCollection<Gene> genes) {
		SAMFileHeader rtrn=new SAMFileHeader();
		SAMSequenceDictionary dict=new SAMSequenceDictionary();
		Iterator<Gene> iter=genes.sortedIterator();
		while(iter.hasNext()){
			Gene gene=iter.next();
			SAMSequenceRecord record=new SAMSequenceRecord(gene.getName(), gene.size());
			dict.addSequence(record);
		}
		rtrn.setSequenceDictionary(dict);
		rtrn.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		return rtrn;
	}

	private Collection<Pair<SAMRecord>> convert(PairedMappedFragment<SAMFragment> read, AnnotationCollection<Gene> genes, SAMFileHeader header) {
		Collection<Pair<SAMRecord>> rtrn=new ArrayList<Pair<SAMRecord>>();
		CloseableIterator<Gene> overlappingGenes=genes.sortedIterator(read, true);
		
		//For each overlapping gene --> convert
		while(overlappingGenes.hasNext()){
			Gene gene=overlappingGenes.next();
			if(read.getOrientation().equals(gene.getOrientation())){
				Pair<SAMRecord> newRead=convert(read, gene, header);
				if(newRead!=null){rtrn.add(newRead);}
			}
		}
		return rtrn;
	}

	private Pair<SAMRecord> convert(PairedMappedFragment<SAMFragment> read, Gene gene, SAMFileHeader header) {
		Annotation newInterval1=gene.convertToFeatureSpace(read.getRead1());
		Annotation newInterval2=gene.convertToFeatureSpace(read.getRead2());
		
		if(newInterval1!=null && newInterval2!=null){
			
			if(newInterval1.size()!=read.getRead1().getSamRecord().getReadString().length() || newInterval2.size()!=read.getRead2().getSamRecord().getReadString().length()){
				System.err.println(read.getRead1().getSamRecord().getReadString()+" "+newInterval1.size()+" "+read.getRead1().getSamRecord().getReadLength()+" "+read.getRead1().getCigarString()+" "+read.getRead1().getSamRecord().toString());
				System.err.println(read.getRead2().getSamRecord().getReadString()+" "+newInterval2.size()+" "+read.getRead2().getSamRecord().getReadLength()+" "+read.getRead2().getCigarString()+" "+read.getRead2().getSamRecord().toString());
				
				return null;
			}
			
			int insertSize=newInterval1.getReferenceEndPosition()-newInterval2.getReferenceStartPosition();
			
			int start1=newInterval1.getReferenceStartPosition()+1;
			int start2=newInterval2.getReferenceStartPosition()+1;
			String sequence1=read.getRead1().getSamRecord().getReadString();
			String sequence2=read.getRead2().getSamRecord().getReadString();
			if(gene.getOrientation().equals(Strand.NEGATIVE)){
				sequence1=Sequence.reverseComplement(sequence1);
				sequence2=Sequence.reverseComplement(sequence2);
				start1++;
				start2++;
			}
			
			SAMRecord record1=new SAMRecord(header);
			record1.setAlignmentStart(start1);
			record1.setCigarString(newInterval1.getCigarString());
			record1.setReferenceName(gene.getName());
			record1.setReadName(read.getName());
			record1.setReadNegativeStrandFlag(read.getRead1().getSamRecord().getReadNegativeStrandFlag());
			record1.setReadString(sequence1);
			record1.setMappingQuality(read.getRead1().getSamRecord().getMappingQuality());
			record1.setBaseQualityString(read.getRead1().getSamRecord().getBaseQualityString());
			record1.setMateReferenceName(newInterval2.getReferenceName());
			record1.setMateAlignmentStart(start2);
			record1.setMateUnmappedFlag(false);
			record1.setInferredInsertSize(insertSize);
			
			
			SAMRecord record2=new SAMRecord(header);
			
			record2.setReferenceName(newInterval2.getReferenceName());
			record2.setAlignmentStart(start2);
			record2.setReadName(read.getName());
			record2.setMateReferenceName(newInterval1.getReferenceName());
			record2.setMateAlignmentStart(start1);
			record2.setMateUnmappedFlag(false);
			record2.setReadString(sequence2);
			record2.setBaseQualityString(read.getRead2().getSamRecord().getBaseQualityString());
			record2.setCigarString(newInterval2.getCigarString());
			record2.setInferredInsertSize(insertSize);
			
			record1.setFirstOfPairFlag(false);
			record2.setFirstOfPairFlag(true);
			record1.setSecondOfPairFlag(true);
			record2.setSecondOfPairFlag(false);
			
			
			record1.setReadNegativeStrandFlag(true);
			record1.setMateNegativeStrandFlag(false);
			record2.setReadNegativeStrandFlag(false);
			record2.setMateNegativeStrandFlag(true);
			
			
			record1.setProperPairFlag(true);
			record2.setProperPairFlag(true);
			
			record1.setReadPairedFlag(true);
			record2.setReadPairedFlag(true);
			
			record1.setMappingQuality(255);
			record2.setMappingQuality(255);
			
			Pair<SAMRecord> newpair=new Pair<SAMRecord>(record1, record2); 
			return newpair;
			
		}
		
		/*if(newInterval!=null){
			if(newInterval.size()!=read.getSamRecord().getReadString().length()){
				System.err.println(read.getSamRecord().getReadString()+" "+newInterval.size()+" "+read.getSamRecord().getReadLength());
				return null;
			}
			
			String sequence=read.getSamRecord().getReadString();
			if(gene.getOrientation().equals(Strand.NEGATIVE)){
				sequence=Sequence.reverseComplement(sequence);
			}
			
			SAMRecord record=newInterval.getSamRecord(header);
			record.setAlignmentStart(newInterval.getReferenceStartPosition()+2);
			record.setCigarString(newInterval.getCigarString());
			record.setReferenceName(gene.getName());
			record.setReadName(read.getName());
			record.setReadNegativeStrandFlag(read.getSamRecord().getReadNegativeStrandFlag());
			record.setReadString(sequence);
			record.setMappingQuality(read.getSamRecord().getMappingQuality());
			record.setBaseQualityString(read.getSamRecord().getBaseQualityString());
			
			
			
			record.setReadPairedFlag(false);
			record.setReadPairedFlag(read.getSamRecord().getReadPairedFlag());
			record.setFirstOfPairFlag(read.getSamRecord().getFirstOfPairFlag());
			record.setSecondOfPairFlag(read.getSamRecord().getSecondOfPairFlag());
			
			record.set
			
			//System.err.println(record.toString());
			
			SAMFragment fragment= new SAMFragment(record);
			//System.err.println(fragment.toBED());
			return fragment;
		}*/
		return null;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File bam=new File(args[0]);
			AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile(args[1]);
			String save=args[2];
			new ConvertReadsToTranscriptome(bam, genes, new File(save));
		} 
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=bam file \n args[1]=BED file \n args[2]=save";
	
}
