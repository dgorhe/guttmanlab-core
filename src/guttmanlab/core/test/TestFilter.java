package guttmanlab.core.test;

import java.io.File;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.predicate.UniqueMapperFilter;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.util.CloseableIterator;

public class TestFilter {

	public static void main(String[] args){
		File sampleBamFile=new File("/Users/mguttman/Desktop/aligned.human.allreads.sorted.bam");
		BAMSingleReadCollection bam=new BAMSingleReadCollection(sampleBamFile);
		bam.addFilter(new UniqueMapperFilter());

		CloseableIterator<SAMFragment> reads=bam.sortedIterator();
		while(reads.hasNext()){
			SAMFragment read=reads.next();
			if(read.getMappingQuality()!=255){
				System.err.println(read.getName()+" "+read.getMappingQuality());
			}
		}
	}
	
}
