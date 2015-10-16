package guttmanlab.core.test;

import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;

import java.io.File;
import java.io.PrintWriter;
import java.util.concurrent.TimeUnit;

import net.sf.samtools.util.CloseableIterator;

public class ReadWriteTest {
	
	//there is a bug in file I/O resulting in a lost character in the name of the reference.
	//program is printing "chr9" instead of "chr19"
	
	public static void main(String args[]) throws Exception{
		final long startTime = System.currentTimeMillis();

		String samplefile = "/Users/cburghard/Downloads/chr19.clean.sorted.bam";//args[0];
		String outfile = "/Users/cburghard/Downloads/readwrite.txt";//args[3];

		//1. read in bamfile
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(new File(samplefile));
		PrintWriter writer = new PrintWriter(new File(outfile));
		
		//2. write bamfile copy
		bamPair.writeToBAM(outfile);
		
		//3. print records from original
		CloseableIterator<PairedMappedFragment<SAMFragment>> bp = bamPair.sortedIterator();
		int count = 0;
		while(bp.hasNext())
		{
			count++;
			if(count > 5) break;
			System.out.println(bp.next().toBED());
		}
		bp.close();
		
		System.out.println("*****************************************");
		
		//4. read and print records from copy
		BAMPairedFragmentCollection bamPairWritten = new BAMPairedFragmentCollection(new File(outfile));
		CloseableIterator<PairedMappedFragment<SAMFragment>> bpw = bamPairWritten.sortedIterator();
		count = 0;
		while(bpw.hasNext())
		{
			count++;
			if (count > 5) break;
			System.out.println(bpw.next().toBED());
		}
		bpw.close();
		writer.close();
		System.out.println("Ran for "+TimeUnit.MINUTES.convert((System.currentTimeMillis()-startTime), TimeUnit.NANOSECONDS)+" min.");
	}
	
}
