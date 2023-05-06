package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.Score;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import guttmanlab.core.datastructures.IntervalTree;
import net.sf.samtools.util.CloseableIterator;

/**
 * 
 * @author mitch
 * A datastructure representing a genome interaction matrix
 * Store a genome in columns and a genome in rows
 * At defined resolution
 */
public class GenomeInteractionMatrix {

	//private BarcodingData data;
	private CoordinateSpace genome;
	private BAMSingleReadCollection bam;
	
	public GenomeInteractionMatrix(File bamFile, CoordinateSpace genome){
		this.bam=new BAMSingleReadCollection(bamFile);
		this.genome=genome;
	}
	
	private Set<String> getBarcodes(BAMSingleReadCollection bam2, Annotation region2) {
		Set<String> rtrn=new TreeSet<String>();
		CloseableIterator<SAMFragment> iter=bam2.sortedIterator(region2, true);
		
		while(iter.hasNext()){
			SAMFragment read=iter.next();
			String name=read.getName();
			String barcode=name.split("::")[1];
			rtrn.add(barcode);
		}
		return rtrn;
	}

	/*public void getTable(int resolution, String save) throws IOException{
		Collection<SingleInterval> regions1=getRegions(resolution);
		writeTable(regions1, regions1, save);
	}*/
	
	/*public void getTable(String chr1, String chr2, int resolution, String save) throws IOException{
		Collection<SingleInterval> regions1 = getRegions(chr1, resolution);
		Collection<SingleInterval> regions2 = getRegions(chr2, resolution);
		writeTable(regions1, regions2, save);
	}*/
	
	/*private void writeTable(Collection<SingleInterval> regions1, Collection<SingleInterval> regions2, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		//writer.write("Name");
		for(SingleInterval interval2:regions2){
			writer.write(interval2.toUCSC()+"\t");
		}
		writer.write("\n");
		
		int counter=0;
		for(SingleInterval interval1: regions1){
			//writer.write(interval1.toUCSC());
			for(SingleInterval interval2: regions2){
				InteractionScore score=getInteractionScore(interval1, interval2);
				
				writer.write(score.getScore()+"\t");
				counter++;
				if(counter % 10000 ==0){System.err.println(counter+" "+interval1.toUCSC());}
			}
			writer.write("\n");
		}
		writer.close();
	}*/

	private Collection<SingleInterval> getRegions(String chr1, int resolution) {
		int chr1Length=genome.getRefSeqLengths().get(chr1);
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		for(int i=0; i<chr1Length; i+=resolution){
			SingleInterval interval=new SingleInterval(chr1, i, i+resolution);
			System.err.println(interval.toUCSC());
			rtrn.add(interval);
		}
		return rtrn;
	}
	
	private Collection<SingleInterval> getRegions(int resolution) {
		Collection<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		for(String chr: genome.getRefSeqLengths().keySet()){
			rtrn.addAll(getRegions(chr, resolution));
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File file=new File(args[0]);
			CoordinateSpace genome=new CoordinateSpace(args[1]);
			
			GenomeInteractionMatrix matrix=new GenomeInteractionMatrix(file, genome);
			
			/*Annotation region1=new SingleInterval("chr1", 10000000,11000000);
			Annotation region2=new SingleInterval("chr1", 12000000,13000000);
			
			System.err.println(matrix.getInteractionScore(region1, region2).getSharedBarcodes());
			*/
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BAM file \n args[1]=sizes \n args[2]=save";
}
