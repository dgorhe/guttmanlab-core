package guttmanlab.core.examples;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import net.sf.samtools.util.CloseableIterator;

public class GencodeExpression {

	public static void main(String[] args) throws IOException,Exception {
		String dir = "/storage/Users/cburghard/Projects/071515_TSS/";
		BEDFileIO io =  new BEDFileIO(dir+"gencode_nonzero.bed");
		BAMSingleReadCollection reads = new BAMSingleReadCollection(new File(dir+"cage.bam"));
		String FeatureFile = dir+"nonexpressed.bed";
		CloseableIterator<? extends Annotation> features = io.loadFromFile(FeatureFile).sortedIterator();
		PrintWriter writer = new PrintWriter("/Users/cburghard/Downloads/geneLengthsNew.txt");
		int count = 0;
		
		while(features.hasNext())
		{
			count++;
			Annotation feature = features.next();
			
			//if(count % 1000 == 0)
			System.out.println(feature.toBED()+"\t"+reads.numOverlappers(feature,false));
			
			if(reads.numOverlappers(feature, false) > 0)
				writer.print(feature.toBED()+"\n");
		}
		
		features.close();
		writer.close();
	}
}
