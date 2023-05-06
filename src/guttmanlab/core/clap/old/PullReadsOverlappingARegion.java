package guttmanlab.core.clap.old;

import java.io.File;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

public class PullReadsOverlappingARegion {

	public PullReadsOverlappingARegion(File bam, Annotation region, File save){
		BAMPairedFragmentCollection bam1=new BAMPairedFragmentCollection(bam);
		SAMFileWriter writer=new SAMFileWriterFactory().makeBAMWriter(bam1.getFileHeader(), false, save);
				
		CloseableIterator<PairedMappedFragment<SAMFragment>> iter=bam1.sortedIterator(region, false);
		while(iter.hasNext()){
			PairedMappedFragment<SAMFragment> reads=iter.next();
			SAMFragment read1=reads.getRead1();
			SAMFragment read2=reads.getRead2();
			writer.addAlignment(read1.getSamRecord());
			writer.addAlignment(read2.getSamRecord());
		}
		
		iter.close();
		writer.close();
	}
	
	public static void main(String[] args){
		if(args.length>2){
			File bam=new File(args[0]);
			File save=new File(args[1]);
			String region=args[2];
			String chr=region.split(":")[0];
			int start=new Integer(region.split(":")[1].split("-")[0]);
			int end=new Integer(region.split(":")[1].split("-")[1]);
			new PullReadsOverlappingARegion(bam, new SingleInterval(chr, start, end), save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam \n args[1]=save \n args[2]=region";
	
}
