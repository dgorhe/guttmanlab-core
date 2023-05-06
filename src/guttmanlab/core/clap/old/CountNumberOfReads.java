package guttmanlab.core.clap.old;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.PopulatedWindow;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import net.sf.samtools.util.CloseableIterator;

public class CountNumberOfReads {

	public CountNumberOfReads(File bamFile, Collection<? extends Annotation> genes, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		BAMSingleReadCollection bam1=new BAMSingleReadCollection(bamFile);
		
		
		
		int counter=0;
		for(Annotation gene: genes){
			SingleInterval genomicWindow=new SingleInterval(gene.getReferenceName(), gene.getReferenceStartPosition(), gene.getReferenceEndPosition(), gene.getOrientation());
			CloseableIterator<? extends PopulatedWindow<SAMFragment>> windows=bam1.getPopulatedWindows(genomicWindow, 100, 100);
			while(windows.hasNext()){
				PopulatedWindow<SAMFragment> window=windows.next();
				
				writer.write(window.toUCSC(gene.getOrientation())+"\t"+window.getNumberOfAnnotationsInWindow()+"\n");
				//System.out.println(window.toUCSC(window.getOrientation())+"\t"+window.getNumberOfAnnotationsInWindow());
			}
			windows.close();
			counter++;
			if(counter % 100 ==0){System.err.println(counter);}
		}
		writer.close();
	}
	
	
	
	private int countMouse(BAMSingleReadCollection bam1) {
		CloseableIterator<SAMFragment> iter=bam1.sortedIterator();
		int counter=0;
		while(iter.hasNext()){
			SAMFragment read=iter.next();
			if(read.getReferenceName().contains("mouse")){counter++;}
		}
		iter.close();
		return counter;
	}



	private static Collection<? extends Annotation> getRegions(Map<String, Collection<Gene>> regions, String chr){
		Collection<Gene> list=new TreeSet<Gene>();
		if(chr.equalsIgnoreCase("ALL")){
			for(String c: regions.keySet()){
				list.addAll(regions.get(c));
			}
		}
		else if(chr.equalsIgnoreCase("human")){
			for(String c: regions.keySet()){
				if(c.contains("human")){
					list.addAll(regions.get(c));
				}
			}
		}
		else if(chr.equalsIgnoreCase("mouse")){
			for(String c: regions.keySet()){
				if(c.contains("mouse")){
					list.addAll(regions.get(c));
				}
			}
		}
		else{
			list=regions.get(chr);
		}
		
		return list;
		
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
			File clipBam=new File(args[0]);
			
			String BEDFile=(args[1]);
			
			String save=args[2];
			
			
			Map<String, Collection<Gene>> regions=BEDFileIO.loadRegionsFromFileByChr(BEDFile);
			
			Collection<? extends Annotation> humanGenes=getRegions(regions, "all");
			
			
			
			new CountNumberOfReads(clipBam, humanGenes, save);
			
			
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=clip BAM file \n args[1]=BED file of regions \n args[2]=save";
	
	
}
