package guttmanlab.core.examples;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import net.sf.samtools.util.CloseableIterator;

/*
 * Author: Christina Burghard
 * Converts a BAM file into a Flat file compliant with the GA4GH API reference server implementation
 */
public class BamToGA4GH {
	
	public static void main(String args[]) throws Exception, IOException{
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam sample file", true);
		p.addStringArg("--id", "Sample id",true);
		p.addStringArg("-g", "Bed gene file", false, "/Volumes/storage/Annotations/RefSeq/mm9/RefSeq.bed");
		p.addStringArg("-s", "Chromsome size file", false, "/Volumes/storage/Users/cburghard/Projects/RAP_Pipeline/mm9chrm.bed");
		p.addStringArg("--score","Score (RPKM/FPKM/TPM)",false,"");
		p.parse(args);
		
		String bamFile = p.getStringArg("-b");
		String id = p.getStringArg("--id");
		String FeatureFile = p.getStringArg("-g");
		
		//TODO: support scoring function
		String units = "reads";
		
		//TODO: figure out how isNormalized is supposed to be set
		boolean isNormalized = false;
		
		//read in bam file
		BAMPairedFragmentCollection bam = new BAMPairedFragmentCollection(new File(bamFile));
		
		//read in feature list
		BEDFileIO io =  new BEDFileIO(p.getStringArg("-s"));
		CloseableIterator<? extends Annotation> features = io.loadFromFile(FeatureFile).sortedIterator();
		//TODO: delete this
		String dir = "/Users/cburghard/Documents/workspace/server/ga4gh-example-data/rnaQuant/diffExpr/";
		PrintWriter exprWriter = new PrintWriter(dir+"expression.table");
		
		//assign an id and annotationId
		//TODO: how to assign appropriate ids?
		String annotationId = "v1";
		
		//TODO: calculate distribution summary statistics
		double complexity = 0;
		double exonicFraction = 0;
		double fractionMapped = 0;
		double intergenicFraction = 0;
		double intronicFraction = 0;
		PrintWriter distWriter = new PrintWriter(dir+"dist.table");
		distWriter.println(id+"\t"+complexity+"\t"+exonicFraction+"\t"+fractionMapped+"\t"+intergenicFraction+"\t"+intronicFraction);
		distWriter.close();
		
		//TODO: calculate counts summary statistics
		int multiCount = 0;
		int multiSpliceCount = 0;
		int totalReadCount = bam.getNumAnnotations();
		int uniqueCount = 0;
		int uniqueSpliceCount = 0;		
		PrintWriter countWriter = new PrintWriter(dir+"counts.table");
		countWriter.println(id+"\t"+multiCount+"\t"+multiSpliceCount+"\t"+totalReadCount+"\t"+uniqueCount+"\t"+uniqueSpliceCount);
		countWriter.close();
		
		//TODO: calculate rnaQuant summary statistics
		PrintWriter quantWriter = new PrintWriter(dir+"rnaseq.table");
		String annotations = "refSeq";
		String description = "a test";
		String name = "mm9";
		String readGroupId = "sample10";
		quantWriter.println(id+"\t"+annotations+"\t"+description+"\t"+name+"\t"+readGroupId);
		quantWriter.close();

		//for each feature, calculate expression score
		Annotation region;	
		String featureGroupId = null;
		int rawReadCount = 0;
		double exprVal = 0;
		double score = 0;
		
		while(features.hasNext())
		{
			region = features.next();
			featureGroupId = region.getName();  //TODO: how to handle lists containing multiple isoforms?
			rawReadCount = bam.numOverlappers(region, true);
			//TODO: correctly calculate expression, score
			exprVal = 0;
			score = 0;
			
			exprWriter.println(id+
				"\t"+annotationId+
				"\t"+exprVal+
				"\t"+featureGroupId+
				"\t"+isNormalized+
				"\t"+rawReadCount+
				"\t"+score+
				"\t"+units);
		}
		
		exprWriter.close();
	}
}
