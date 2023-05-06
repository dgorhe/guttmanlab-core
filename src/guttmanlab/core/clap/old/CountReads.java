package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.util.CloseableIterator;

public class CountReads {
	
	public CountReads(BAMSingleReadCollection reads, Collection<SingleInterval> regions, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: regions){
			int count=reads.numOverlappers(region, true);
			int countStranded=0;
			//TODO Count iff it starts within the region
			CloseableIterator<SAMFragment> iter=reads.sortedIterator(region, true);
			while(iter.hasNext()){
				SAMFragment read=iter.next();
				if(read.getOrientation().equals(Strand.NEGATIVE)){countStranded++;}
			}
			iter.close();
			writer.write(region.toUCSC()+"\t"+count+"\t"+countStranded+"\n");
		}
		
		writer.close();
	}
	
	
	public static void main (String[] args) throws IOException{
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(new File(args[0]));
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		String save=args[2];
		new CountReads(bam1, regions, save);
	}
	
}
