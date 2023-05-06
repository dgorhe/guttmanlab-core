package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.datastructures.Pair;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class GetReadsWithinSNPPosition {

	public GetReadsWithinSNPPosition(File bamFile, VCFFileReader snps, String fileName, String sp1, String sp2){
		SAMFileReader reader=new SAMFileReader(bamFile);
		
		SAMRecordIterator reads=reader.iterator();
		
		Collection<String> list=new TreeSet<String>();
		
		int counter=0;
		int specific=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			boolean containsSNP=overlapsSNP(record, snps, sp1, sp2);
			if(containsSNP){
				list.add(record.getReadName());
				specific++;
			}
			
			counter++;
			if(counter%10000 ==0){System.err.println(counter+" "+specific);}
		}
		
		
		reads.close();
		
		write(reader, reader.getFileHeader(), list, fileName);
		
		
		reader.close();
		reads.close();
		snps.close();
	}

	private boolean overlapsSNP(SAMRecord record, VCFFileReader snps, String sp1, String sp2) {
		CloseableIterator<VariantContext> iter=snps.query(record.getReferenceName().replace("chr", ""), record.getAlignmentStart(), record.getAlignmentEnd());
		while(iter.hasNext()){
			VariantContext variant=iter.next();
			boolean isHet=isHetero(variant, sp1, sp2);
			
			if(isHet){
				iter.close();
				return true;
			}
		}
		
		iter.close();
		return false;
	}
	
	private boolean isHetero(VariantContext variant, String sp1, String sp2) {
		return !variant.getGenotype(sp1).sameGenotype(variant.getGenotype(sp2));
	}

	private void write(SAMFileReader reader, SAMFileHeader fileHeader, Collection<String> list, String fileName) {
		SAMFileWriter writer1=new SAMFileWriterFactory().setCreateIndex(true).makeSAMOrBAMWriter(fileHeader, false, new File(fileName));
		
		SAMRecordIterator reads=reader.iterator();
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			if(list.contains(record.getReadName())){
				writer1.addAlignment(record);
			}
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		writer1.close();
		reads.close();

	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			
		File file=new File(args[0]);
		VCFFileReader vcf=new VCFFileReader(new File(args[1]));
		String fileName=args[2];
		String species1=args[3];
		String species2=args[4];
		new GetReadsWithinSNPPosition(file, vcf, fileName, species1, species2);
		}
		else{System.err.println(usage);}
		
		
	}
	
	static String usage=" args[0]=bam file \n args[1]=VCF file \n args[2]=save file name \n args[3]=species name 1 \n args[4]=species name 2";

	
}
