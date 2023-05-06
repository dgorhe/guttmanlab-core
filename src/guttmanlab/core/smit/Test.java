package guttmanlab.core.smit;

import java.io.File;

import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.util.CloseableIterator;

public class Test {

	public static void main(String[] args){
		/*BAMSingleReadCollection single=new BAMSingleReadCollection(new File(args[0]));
		
		CloseableIterator<SAMFragment> reads2=single.sortedIterator();
		while(reads2.hasNext()){
			SAMFragment read=reads2.next();
			if(read.getName().contains("M05340:38:000000000-BL6FL:1:1109:7068:8660")){
				System.err.println("single data has "+read.getName()+"\n"+read.toString());
			}
		}
		reads2.close();*/
		
		
		BAMPairedFragmentCollection paired=new BAMPairedFragmentCollection(new File(args[0]));
		CloseableIterator<PairedMappedFragment<SAMFragment>> reads=paired.sortedIterator();
		while(reads.hasNext()){
			PairedMappedFragment<SAMFragment> fragment=reads.next();
			if(fragment.getName().contains("M05340:38:000000000-BL6FL:1:1109:7068:8660")){
				System.err.println("paired data has "+fragment.getName());
			}
		}
		reads.close();
		
		
	}
	
}
