package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.BlockedAnnotation;
import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;

public class Triples {

	public Triples(Map<String, IntervalTree<File>> barcodeFileTree, Cluster triple, int resolution, String save) throws IOException{
		SingleInterval region=genomeSpan(triple);
		Iterator<File> iter=barcodeFileTree.get(region.getReferenceName()).overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
		BarcodingData data=new BarcodingData(iter);
		
		int observed=data.quantify(triple);
		
		System.err.println("observed "+triple.toStringNoName()+" "+observed);
		
		//TODO Fix i,j --> vary k at fixed resolution
		
		ArrayList<SingleInterval> bins=bins(region, resolution);
		Iterator<SingleInterval> iter1=triple.getAllIntervals().iterator();
		SingleInterval interval1=iter1.next();
		SingleInterval interval2=iter1.next();
		SingleInterval interval3=iter1.next();
		
		Map<Cluster, Integer> scores=new TreeMap<Cluster, Integer>();
		
		for(SingleInterval bin: bins){
			ArrayList<SingleInterval> bins2=bins(interval2, resolution);
			ArrayList<SingleInterval> bins3=bins(interval3, resolution);
			for(SingleInterval i2: bins2){
				for(SingleInterval i3: bins3){
					Cluster random=new Cluster("r");
					random.addRead(bin);
					random.addRead(i2);
					random.addRead(i3);
					int score=data.quantify(random);
					scores.put(random, score);
					System.err.println(random.toStringNoName()+" "+score);
				}
			}
		}
		
		for(SingleInterval bin: bins){
			ArrayList<SingleInterval> bins1=bins(interval1, resolution);
			ArrayList<SingleInterval> bins3=bins(interval3, resolution);
			
			for(SingleInterval i1: bins1){
				for(SingleInterval i3: bins3){
					Cluster random=new Cluster("r");
					random.addRead(i1);
					random.addRead(bin);
					random.addRead(i3);
					int score=data.quantify(random);
					scores.put(random, score);
					System.err.println(random.toStringNoName()+" "+score);
				}
			}
			
		}
		
		for(SingleInterval bin: bins){
			ArrayList<SingleInterval> bins1=bins(interval1, resolution);
			ArrayList<SingleInterval> bins2=bins(interval2, resolution);
			
			for(SingleInterval i1: bins1){
				for(SingleInterval i2: bins2){
					Cluster random=new Cluster("r");
					random.addRead(i1);
					random.addRead(i2);
					random.addRead(bin);
					int score=data.quantify(random);
					scores.put(random, score);
					System.err.println(random.toStringNoName()+" "+score);
				}
			}
			
			
		}
		
		
		
		/*Iterator<SingleInterval> iter1=triple.getAllIntervals().iterator();
		ArrayList<SingleInterval> bins1=bins(region, iter1.next().getLength());
		ArrayList<SingleInterval> bins2=bins(region, iter1.next().getLength());
		ArrayList<SingleInterval> bins3=bins(region, iter1.next().getLength());*/
		
		
		
		
		/*int counter=0;
		//Enumerate ALL triples in region
		for(SingleInterval i1: bins1){
			for(SingleInterval i2: bins2){
				for(SingleInterval i3: bins3){
					if(!i1.overlaps(i2) && !i1.overlaps(i3) && !i2.overlaps(i3)){
						Cluster random=new Cluster("r");
						random.addRead(i1);
						random.addRead(i2);
						random.addRead(i3);
						int score=data.quantify(random);
						scores.put(random, score);
						counter++;
						if(counter%1 ==0){System.err.println(counter+" "+i1.toUCSC()+" "+i2.toUCSC()+" "+i3.toUCSC()+" "+score);}
					}
				}
			}
			
		}*/
		
		write(save, scores);
		
	}

	private ArrayList<SingleInterval> bins(SingleInterval region, int resolution) {
		Iterator<DerivedAnnotation<? extends Annotation>> iter=region.getWindows(resolution, resolution).sortedIterator();
		ArrayList<SingleInterval> rtrn=new ArrayList<SingleInterval>();
		
		while(iter.hasNext()){
			DerivedAnnotation<? extends Annotation> a=iter.next();
			SingleInterval SI=new SingleInterval(a.getReferenceName(), a.getReferenceStartPosition(), a.getReferenceEndPosition());
			rtrn.add(SI);
		}
		
		return rtrn;
	}

	private SingleInterval genomeSpan(Cluster triple) {
		SingleInterval interval=triple.getAllIntervals().iterator().next();
		String chr=interval.getReferenceName();
		int start=interval.getReferenceStartPosition();
		int end=interval.getReferenceEndPosition();
		
		for(SingleInterval i: triple.getAllIntervals()){
			chr=i.getReferenceName();
			start=Math.min(start, i.getReferenceStartPosition());
			end=Math.max(end, i.getReferenceEndPosition());
		}
		return new SingleInterval(chr, start, end);
	}

	private void write(String save, Map<Cluster, Integer> scores) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(Cluster c: scores.keySet()){
			int count=scores.get(c);
			Annotation blocks=new BlockedAnnotation(c.getAllIntervals(), "s="+count);
			if(blocks.getNumberOfBlocks()==3){
				writer.write(blocks.toBED()+"\n");
			}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		Map<String, IntervalTree<File>> barcodeFileTree=PermuteKmers.makeFileTree(files);
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		int resolution=new Integer(args[2]);
		String save=args[3];
		Cluster c=new Cluster("triple");
		c.addReads(regions);
		new Triples(barcodeFileTree, c, resolution, save);
	}
	
}
