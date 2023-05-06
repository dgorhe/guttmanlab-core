package guttmanlab.core.xist;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class CountXToAUnique {

	private static void getUniqueReads(File bam, String save) throws IOException {
		SAMFileReader reader=new SAMFileReader(bam);
		SAMRecordIterator reads=reader.iterator();
		
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			//if(!record.getFirstOfPairFlag()) {
				SingleInterval interval=new SingleInterval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentStart()+1);
				rtrn.add(interval);
			//}
			
			counter++;
			if(counter%1000000 ==0){System.err.println(counter);}
		}
		
		reader.close();
		reads.close();
		
		BEDFileIO.writeBED(rtrn, save);
		
		double Xcount=0;
		double totalCount=0;
		for(SingleInterval r: rtrn) {
			if(r.getReferenceName().equals("chrX")) {Xcount++;}
			totalCount++;
		}
		double XA=Xcount/totalCount;
		
		System.err.println(Xcount+" "+totalCount+" "+XA);
		
	}
	
	public static void main(String[] args) throws IOException {
		File bam=new File(args[0]);
		String save=args[1];
		getUniqueReads(bam, save);
	}
}
