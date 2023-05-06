package guttmanlab.core.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.math.Statistics;
import net.sf.samtools.util.CloseableIterator;

public class ComputeEnrichmentPerWindow {

	
	public static void main(String[] args) throws IOException{
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(new File(args[0]));
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		AnnotationCollection<Gene> genes=BEDFileIO.loadFromFile(args[2]);
		String save=args[3];
		int windowLength=new Integer(args[4]);
		
		FileWriter writer=new FileWriter(save);
		
		
		for(SingleInterval region: regions){
			//for each region get all 5Kb windows
			CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows=bam1.getPopulatedWindows(region, windowLength, windowLength);
			
			List<Double> scores=new ArrayList<Double>();
			int counter=0;
			while(windows.hasNext()){
				PopulatedWindow<SAMFragment> window=windows.next();
				double count=window.getNumberOfAnnotationsInWindow();
				double normPerBase=count/(double)windowLength;
				scores.add(count);
				counter++;
				if(counter%10 ==0){System.err.println(counter+" "+window.toUCSC()+" "+count+" "+normPerBase+" "+getDensity(window, genes));}
				writer.write(window.toUCSC()+"\t"+count+"\t"+getDensity(window, genes)+"\n");
			}
			windows.close();
			
			Collections.sort(scores);
			
			double average=Statistics.mean(scores);
			double median=Statistics.quantile(scores, 0.5);
			double p90=Statistics.quantile(scores, 0.9);
			double geneDensity=getDensity(region, genes);
			
			System.out.println(region.toUCSC()+"\t"+average+"\t"+median+"\t"+p90+"\t"+geneDensity);
		}
		
		writer.close();
		
	}

	private static double getDensity(Annotation region, AnnotationCollection<Gene> genes) {
		/*double regionLength=(double)region.getReferenceEndPosition()-region.getReferenceStartPosition();
		return (double)genes.numOverlappers(region, false)/regionLength;*/
		return genes.numOverlappers(region, false);
	}
	
}
