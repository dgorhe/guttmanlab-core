package guttmanlab.core.test;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PairedMappedFragment;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.util.CommandLineParser;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;

import org.apache.commons.lang.time.StopWatch;

import jsc.distributions.Binomial;
import net.sf.samtools.util.CloseableIterator;
import guttmanlab.core.math.ScanStat;

public class SimplePeakCaller {
	public static void main(String[] args) throws IOException
	{
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam sample file", true);
		p.addStringArg("-i", "Bam input file", true);
		p.addStringArg("-o", "Output prefix",true);
		p.addStringArg("-g", "Bed gene file", false, "/storage/Annotations/RefSeq/mm9/RefSeq.bed");
		p.addStringArg("-s", "Chromsome size file", false, "/storage/Users/cburghard/Projects/RAP_Pipeline/mm9chrm.bed");
		p.addIntArg("-ws", "Window size", false, 5000);
		p.addIntArg("-step", "Window step size", false, 500);
		p.addDoubleArg("-p", "P value cutoff", false, .005);
		p.parse(args);
		
		String SampleFile = p.getStringArg("-b");
		String ControlFile = p.getStringArg("-i");
		String OutputFile = p.getStringArg("-o");
		String FeatureFile = p.getStringArg("-g");
		int winSize = p.getIntArg("-ws");
		int stepSize = p.getIntArg("-step");
		double maxP = p.getDoubleArg("-p");
		
		PrintWriter writer = new PrintWriter(OutputFile, "UTF-8");
		BEDFileIO io =  new BEDFileIO(p.getStringArg("-s"));
		CloseableIterator<? extends Annotation> features = io.loadFromFile(FeatureFile).sortedIterator();
		Annotation region = null;
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(new File(SampleFile));
		BAMPairedFragmentCollection bamInput = new BAMPairedFragmentCollection(new File(ControlFile));
		
		while(features.hasNext())
		{
			
			region = features.next();
			System.out.println("Calculating windows for "+region.getName()); //getName());
			
			int controlTotal = bamInput.numOverlappers(region, false);
			int sampleTotal =  bamPair.numOverlappers(region, false);

			CloseableIterator<? extends PopulatedWindow<PairedMappedFragment<SAMFragment>>> sampleWindows = bamPair.getPopulatedWindows(region, winSize,stepSize,true);
			CloseableIterator<? extends PopulatedWindow<PairedMappedFragment<SAMFragment>>> controlWindows = bamInput.getPopulatedWindows(region, winSize,stepSize,true);

			PopulatedWindow<? extends PairedMappedFragment<SAMFragment>> currentSampleWin = null;
			PopulatedWindow<? extends PairedMappedFragment<SAMFragment>> currentControlWin = null;
			
			while(sampleWindows.hasNext() && controlWindows.hasNext())
			{
				currentSampleWin = sampleWindows.next();
				currentControlWin = controlWindows.next();
				int start = currentSampleWin.getReferenceStartPosition();
				
				int sampleCount = currentSampleWin.getNumberOfAnnotationsInWindow();
				int controlCount = currentControlWin.getNumberOfAnnotationsInWindow();
				double pValue = ScanStat.getPValue(controlCount,sampleCount,controlTotal,sampleTotal,winSize,region.size());
				if(pValue <= maxP)
					{
						writer.println(region.getName()+"\t"+start+"\t"+(start+winSize)+"\t"+pValue);
						//System.out.println(region.getName()+"\t"+start+"\t"+(start+winSize)+"\t"+pValue);
					}
			}
		}	

		writer.close();
	}
	
}