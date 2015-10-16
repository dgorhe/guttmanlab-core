package guttmanlab.core.integrationTests;

import java.io.IOException;
import java.io.PrintWriter;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.io.BEDFileIO;
import net.sf.samtools.util.CloseableIterator;

public class refSeqGeneLengths {

	public static void main(String[] args) throws IOException {
		BEDFileIO io =  new BEDFileIO("/Users/cburghard/Downloads/sizes");
		String FeatureFile = "/Users/cburghard/Downloads/RefSeqNew.bed";
		CloseableIterator<? extends Annotation> features = io.loadFromFile(FeatureFile).sortedIterator();
		PrintWriter writer = new PrintWriter("/Users/cburghard/Downloads/geneLengthsNew.txt");
		
		while(features.hasNext())
		{
			Annotation feature = features.next();
			writer.println(feature.getName()+"\t"+feature.size());
		}
		
		features.close();
		writer.close();
	}

}
