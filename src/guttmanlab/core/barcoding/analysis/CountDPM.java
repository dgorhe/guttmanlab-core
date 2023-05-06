package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.fastq.FastqRecord;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

public class CountDPM {

	public static void main(String[] args) throws IOException{
		
		SAMFileReader bam=new SAMFileReader(new File(args[0]));
		Collection<String> odd=parse(args[1]);
		
		SAMRecordIterator recordIter=bam.iterator();
		
		int countDNA=0;
		int countTotal=0;
		while(recordIter.hasNext()){
			SAMRecord record=recordIter.next();
			//System.err.println(record.getCigarString()+" "+record.getReadString());
			String clippedSequence=getClippedSequence(record);
			/*if(record.getCigarString().contains("D")){
				System.err.println(record.getCigarString()+" "+record.getReadString()+" "+clippedSequence);
			}*/
			
			if(overlaps(clippedSequence, odd)){
				countDNA++;
			}
			
			/*String cigar=record.getCigarString();
			if(cigar.contains("S")){
				//System.err.println(record.getFirstOfPairFlag()+" "+cigar);
				if(record.getFirstOfPairFlag()){
					if(record.getReadString().contains("CAAGTCA")){
						countDNA++;
					}
				}
				else if(!record.getFirstOfPairFlag()){
					if(record.getReadString().contains("TGACTT")){
						countDNA++;
					}
				}
				//countDNA++;
			}*/
			countTotal++;
		}
		
		double ratio=(double)countDNA/(double)countTotal;
		System.err.println(countDNA+" "+countTotal+" "+ratio);
		
		bam.close();
		
	}
	
	
	
	private static boolean overlaps(String clippedSequence, Collection<String> odd) {
		if(clippedSequence.isEmpty() || clippedSequence.equals("")){return false;}
		
		if(clippedSequence.length()<8){return false;}
		
		for(String bc: odd){
			if(bc.contains(clippedSequence)){
				System.err.println(bc+" "+clippedSequence);
				return true;
			}
		}
		
		return true;
	}



	private static String getClippedSequence(SAMRecord record) {
		List<CigarElement> cigar=record.getCigar().getCigarElements();
		
		int countSoFar=0;
		for(CigarElement e: cigar){
			//System.err.println(e.getOperator());
			if(e.getOperator().toString().equals("D") || e.getOperator().toString().equals("N")){
				//countSoFar-=e.getLength();
				//DO NOTHING	
			}
			
			else{
				if(e.getOperator().toString().equals("S")){
					int length=e.getLength();
					return record.getReadString().substring(countSoFar, countSoFar+length);
				}
				countSoFar+=e.getLength();
			}
		}
		return "";
	}



	/*public static void main(String[] args) throws IOException{
		FastqReader reader=new FastqReader(new File(args[0]));
		Collection<String> dpm=parse(args[1]);
		String save=args[2];
		
		
		Collection<String> dnaNames=new TreeSet<String>();
	
		
		SAMFileHeader mouseheader=CoordinateSpace.MM9.getSAMFileHeaderOldVersion();
		SAMFileHeader humanheader=CoordinateSpace.HG19.getSAMFileHeaderOldVersion();
		
		FileWriter dnaWriterH=new FileWriter(save+".human_dna.sam");
		FileWriter rnaWriterH=new FileWriter(save+".human_rna.sam");
		
		writeHeader(dnaWriterH, humanheader);
		writeHeader(rnaWriterH, humanheader);
		
		FileWriter dnaWriterM=new FileWriter(save+".mouse_dna.sam");
		FileWriter rnaWriterM=new FileWriter(save+".mouse_rna.sam");
		
		writeHeader(dnaWriterM, mouseheader);
		writeHeader(rnaWriterM, mouseheader);
		
		
		
		int has=0;
		int total=0;
		//FileWriter writerDNA=new FileWriter(new File(save+"_dna.fq"));
		//FileWriter writerRNA=new FileWriter(new File(save+"_rna.fq"));
		Iterator<FastqRecord> iter=reader.iterator();
		while(iter.hasNext()){
			FastqRecord record=iter.next();
			boolean hasDpm=recordHas(record, dpm);
			if(hasDpm){
				dnaNames.add(record.getReadHeader().split(" ")[0]);
				has++;
				//write(record, writerDNA);
			}
			else{
				//write(record, writerRNA);
			}
			total++;
		}
		double ratio=(double)has/(double)total;
		System.err.println(has+" "+total+" "+ratio);
		reader.close();
		
		
		SAMFileReader bam=new SAMFileReader(new File(args[3]));
		
		SAMRecordIterator recordIter=bam.iterator();
		
		while(recordIter.hasNext()){
			SAMRecord record=recordIter.next();
			String name=record.getReadName().split("::")[0];
			if(dnaNames.contains(name)){
				write(dnaWriterH, dnaWriterM, record);
			}
			else{
				write(rnaWriterH, rnaWriterM, record);
			}
		}
		
		bam.close();
		dnaWriterH.close();
		dnaWriterM.close();
		rnaWriterH.close();
		rnaWriterM.close();
	}*/
	
	private static void writeHeader(FileWriter writer, SAMFileHeader mouseheader) throws IOException {
		List<SAMSequenceRecord> list=mouseheader.getSequenceDictionary().getSequences();
		for(SAMSequenceRecord record: list){
			writer.write("@SQ\tSN:"+record.getSequenceName()+"\tLN:"+record.getSequenceLength()+"\n");
		}
		
		List<String> comments=mouseheader.getComments();
		for(String co: comments){
			writer.write(co+"\n");
			System.err.println(co);
		}
	
	}

	private static void write(FileWriter dnaWriterH, FileWriter dnaWriterM, SAMRecord record) throws IOException {
		FileWriter writer=dnaWriterM;
		//System.err.println(record.getReadName());
		if(record.getReferenceName().contains("human")){
			writer=dnaWriterH;
		}
		
		
		String chr=record.getReferenceName();
		String newChr=trim(chr);
		
		//String mateChr=record.getMateReferenceName();
		//String newMateChr=trim(mateChr);
		
		/*record.setReferenceName(newChr);
		record.setMateReferenceName(newMateChr);*/
		
		if(newChr.contains("chr")){
			writer.write(record.getSAMString());
		}
	}

	private static String trim(String chr) {
		return chr.replace("_mouse", "").replace("_human", "").trim();
	}

	private static void write(FastqRecord record, FileWriter writer) throws IOException {
		writer.write(record.getReadHeader()+"\n");
		writer.write(record.getReadString()+"\n");
		writer.write("+\n");
		writer.write(record.getBaseQualityString()+"\n");
	}

	private static Collection<String> parse(String file) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		List<String> lines=BEDFileIO.loadLines(file);
		for(String line: lines){
			rtrn.add(line.split("\t")[2]);
		}
		return rtrn;
	}

	private static boolean recordHas(FastqRecord record, Collection<String> dpm) {
		for(String seq: dpm){
			String subsequence=record.getReadString().substring(0, 8);
			if(subsequence.equals(seq)){return true;}
		}
		return false;
	}
	
	
	
}
