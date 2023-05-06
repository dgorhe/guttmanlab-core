package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.util.CloseableIterator;

public class WindowCounts {

	private int windowSize;
	private int counter;
	private int readCount;

	public WindowCounts(File bamFile, Collection<Annotation> regions, int windowSize) throws IOException{
		this.windowSize=windowSize;
		this.counter=0;
		computeEnrichment(bamFile, regions);
	}
	
	
	public WindowCounts(File bamFile, int windowSize){
		this.windowSize=windowSize;
		this.counter=0;
		this.readCount=0;
		computeEnrichment(bamFile);
	}
	
	private void computeEnrichment(File bamFile) {
		BAMSingleReadCollection bam=new BAMSingleReadCollection(bamFile);
		//FileWriter writer=new FileWriter(save);
		
		int index=0;
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows1=bam.getPopulatedWindows(windowSize, windowSize);
		while(windows1.hasNext()){
			PopulatedWindow<SAMFragment> window=windows1.next();
			int count=window.getNumberOfAnnotationsInWindow();
			this.readCount+=count;
			if(count>0){
				this.counter++;
				//writer.write(window.toBED(count)+"\n");
			}
			index++;
			//if(index%100000 ==0){System.err.println(index+" "+window.toUCSC()+" "+count);}
		}
		windows1.close();
		bam.close();
		//writer.close();
		
	}
	
	
	private void computeEnrichment(File bamFile, Collection<Annotation> regions) throws IOException {
		BAMSingleReadCollection bam=new BAMSingleReadCollection(bamFile);
		
		int count=0;
		
		for(Annotation a: regions){
			computeEnrichment(bam, a.getOrientation(), a);
			count++;
			if(count % 1000 ==0){System.err.println(a.getName()+" "+a.toUCSC());}
		}
		
		for(Annotation a: regions){
			Collection<Annotation> introns=a.getIntrons();
			for(Annotation intron: introns){
				computeEnrichment(bam, intron.getOrientation(), intron);
			}
			count++;
			if(count %10000 ==0){System.err.println("introns "+a.getName());}
		}
		
		
	}
	
	private void computeEnrichment(BAMSingleReadCollection bam, Strand strand, Annotation region) throws IOException {
		//Iterate over windows
		CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows1=bam.getPopulatedWindows(region, windowSize, windowSize);
		
		
		while(windows1.hasNext()){
			PopulatedWindow<SAMFragment> window=windows1.next();
			window.setOrientation(strand);
			int count=window.getNumberOfAnnotationsInWindow();
			if(count>0){
				this.counter++;
				//writer.write(window.toUCSC(strand)+"\t"+count+"\n");
			}
		}
		windows1.close();
	}
	
	public int getNumberOfWindows(){return this.counter;}
	
	public int getTotalNumberOfReads() {return readCount;}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File bam=new File(args[0]);
			int windowSize=new Integer(args[1]);
			String save=args[2];
			
			WindowCounts counter=new WindowCounts(bam, windowSize);
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam file \n args[1]=window size \n args[2]=save";

	
}
