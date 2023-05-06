package guttmanlab.core.spidr;

import java.io.File;
import java.util.Iterator;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class FindPolyARuns {

	//Read2 would start with TTTTT
	//Read1 would end with AAAA
	
	public FindPolyARuns(File read1, File read2, String seqMatch) {
		FastqReader r2=new FastqReader(read2);
		FastqReader r1=new FastqReader(read1);
		Iterator<FastqRecord> iter2=r2.iterator();
		Iterator<FastqRecord> iter1=r1.iterator();
		while(iter2.hasNext()) {
			//iter1.hasNext();
			FastqRecord record1=iter1.next();
			FastqRecord record2=iter2.next();
			//System.err.println(record1.getReadHeader()+" "+record2.getReadHeader());
			String seq=record2.getReadString();
			if(seq.endsWith(seqMatch)) {
				System.out.println(">"+record1.getReadHeader()+"\n"+record1.getReadString());
			}
		}
		r1.close();
		r2.close();
	}
	
	public static void main(String[] args) {
		File r1=new File(args[0]);
		File r2=new File(args[1]);
		String seqMatch=args[2];
		new FindPolyARuns(r1, r2, seqMatch);
	}
}
