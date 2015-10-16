package guttmanlab.core.examples;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;
import guttmanlab.core.util.CommandLineParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.concurrent.TimeUnit;

public class SimpleNormalize {
	
	//can we read a file on the mounted drive?
	public static void main(String args[]) throws Exception{
		final long startTime = System.currentTimeMillis();

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam sample file", true);
		p.addStringArg("-i", "Bam input file", true);
		p.addStringArg("-o", "Output file",true);
		p.addStringArg("-s", "Chromosome size file", true);
		p.addIntArg("-ws", "window size", false,10000);
		p.parse(args);
		
		//read test
		int win_size = p.getIntArg("-ws");
		int step_size = p.getIntArg("-ws");
		String sizefile = p.getStringArg("-s");
		String samplefile = p.getStringArg("-b");
		String inputfile = p.getStringArg("-i");
		String outfile = p.getStringArg("-o");
		
		BAMPairedFragmentCollection bamPair = new BAMPairedFragmentCollection(new File(samplefile));
		BAMPairedFragmentCollection bamInput = new BAMPairedFragmentCollection(new File(inputfile));
		HashMap<String,Integer> chrs = ChrMapFromSizeFile(new File(sizefile));
		PrintWriter writer = new PrintWriter(new File(outfile));
		
		int st = bamPair.getNumAnnotations();
		int it = bamInput.getNumAnnotations();
		
		for ( String chr : chrs.keySet())
		{
			for (int i=0; i < chrs.get(chr)-(win_size+step_size); i+=step_size)
			{
				SingleInterval region = new SingleInterval(chr,i,i+win_size,Strand.BOTH);
				int sc = (bamPair.numOverlappers(region, false));
				int ic = (bamInput.numOverlappers(region, false));
				double norm_val = normalize(sc,ic,st,it);
				writer.println(chr+"\t"+i+"\t"+(i+win_size)+"\t"+norm_val);
			}
		}
		writer.close();
		System.out.println("Ran for "+TimeUnit.MINUTES.convert((System.currentTimeMillis()-startTime), TimeUnit.MILLISECONDS)+" min.");
	}
	
	private static HashMap<String,Integer> ChrMapFromSizeFile(File f) throws FileNotFoundException, IOException
	{
		HashMap<String,Integer> chrs = new HashMap<String,Integer>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		
		while((line = br.readLine()) != null)
		{
			String[] parts = line.split("\t");
			String chr = parts[0];
			Integer size = Integer.valueOf(parts[2]);
			chrs.put(chr, size-1);
		}
		br.close();
		return chrs;
	}
	
	public static double normalize(int sampleCount, int inputCount, int sampleTotal, int inputTotal)
	{
		double s = (double)sampleCount/(double)sampleTotal;
		double i = (double)inputCount/(double)inputTotal;
		//System.out.println(s+"\t"+i+"\t"+s/i);
		if(i!=0)
			return s/i;  //is this the correct behavior?
		return 0;	
	}
}
