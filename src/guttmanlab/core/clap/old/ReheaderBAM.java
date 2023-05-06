package guttmanlab.core.clap.old;

import java.io.File;
import java.util.List;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.util.CloseableIterator;

public class ReheaderBAM {

	public ReheaderBAM(File bamWithHeader, File newBamFile, String save){
		BAMSingleReadCollection oldBam=new BAMSingleReadCollection(bamWithHeader);
		BAMSingleReadCollection newBam=new BAMSingleReadCollection(newBamFile);
		
		List<String> comments=oldBam.getFileHeader().getComments();
		
		SAMFileHeader header=newBam.getFileHeader();
		for(String comment: comments){
			header.addComment(comment);
		}
		
		BAMFileWriter writer=new BAMFileWriter(new File(save));
		writer.setHeader(header);
		
		CloseableIterator<SAMFragment> reads=newBam.sortedIterator();
		while(reads.hasNext()){
			SAMFragment read=reads.next();
			writer.addAlignment(read.getSamRecord());
		}
		
		reads.close();
		writer.close();
	}
	
	
	public static void main(String[] args){
		if(args.length>2){
			new ReheaderBAM(new File(args[0]), new File(args[1]), args[2]);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam with header \n args[1]=bam without header \n args[2]=new headered file";
	
}
