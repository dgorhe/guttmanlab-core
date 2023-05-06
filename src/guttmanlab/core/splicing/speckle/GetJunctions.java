package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;

public class GetJunctions {

	
	public static void main(String[] args) throws IOException {
		Collection<Gene> transcripts=GTFToJunctions.getTranscripts(new File(args[0]));
		
		Collection<Junction> junctionSet=new TreeSet<Junction>();
		
		for(Gene gene: transcripts) {
			Collection<Gene> junctions=gene.getJunctions();
			for(Gene j: junctions) {
				j.setName(gene.getName());
				int tss=gene.getReferenceStartPosition();
				Junction junc=new Junction(j, tss);
				junctionSet.add(junc);
			}	
		}
		
		write(junctionSet);
		
	}

	private static void write(Collection<Junction> junctionSet) {
		for(Junction j: junctionSet) {
			System.out.println(j.intron.toShortBED()+"\t"+"tss="+j.tss.getReferenceStartPosition());
		}
		
	}

	private static Collection<Gene> filterSingles(Map<String, IntervalTree<Gene>> junctionTree, Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(Gene g: genes) {
			if(g.getNumberOfBlocks()==1) {
				Collection<Gene> overlappers=getOverlappingSingles(g, junctionTree);
				rtrn.addAll(overlappers);
			}
		}
		
		return rtrn;
	}

	private static Collection<Gene> getOverlappingSingles(Gene g, Map<String, IntervalTree<Gene>> junctionTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		if(junctionTree.containsKey(g.getReferenceName())) {
			IntervalTree<Gene> tree=junctionTree.get(g.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(g.getReferenceStartPosition(), g.getReferenceEndPosition());
			while(iter.hasNext()) {
				Gene o=iter.next();
				if(o.getOrientation().equals(g.getOrientation())) {
						rtrn.add(o);
					}
				
			}
		}
		
		return rtrn;
	}

	private static Collection<Gene> getOverlappers(Gene junc, Map<String, IntervalTree<Gene>> junctionTree) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		if(junctionTree.containsKey(junc.getReferenceName())) {
			SingleInterval j=new SingleInterval(junc.getReferenceName(), junc.getFirstBlock().getReferenceEndPosition(), junc.getLastBlock().getReferenceStartPosition());
			IntervalTree<Gene> tree=junctionTree.get(junc.getReferenceName());
			Iterator<Gene> iter=tree.overlappingValueIterator(j.getReferenceStartPosition(), j.getReferenceEndPosition());
			while(iter.hasNext()) {
				Gene o=iter.next();
				SingleInterval oj=new SingleInterval(o.getReferenceName(), o.getFirstBlock().getReferenceEndPosition(), o.getLastBlock().getReferenceStartPosition());
				if(!o.equals(junc)) {
					if(o.getOrientation().equals(junc.getOrientation())) {
						if(j.overlaps(oj)) {
							rtrn.add(o);
						}
					}
				}
			}
		}
		
		return rtrn;
	}
	
}
