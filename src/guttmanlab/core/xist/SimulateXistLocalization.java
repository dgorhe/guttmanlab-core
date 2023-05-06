package guttmanlab.core.xist;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.simulation.CoordinateSpace;

public class SimulateXistLocalization {

	int numRNACopies=200;
	
	public SimulateXistLocalization(String save, int numPerms, IntervalTree<Gene> genes) throws IOException {
		SingleInterval chrX=new SingleInterval("chrX", 0, CoordinateSpace.MM9.getRefSizes().get("chrX"));
		
		List<String> geneNames=getRow(genes);
		List<String> columns=getColumns(numPerms);
		MatrixWithHeaders genesByCell=new MatrixWithHeaders(geneNames, columns);
		
		for(int i=0; i<numPerms; i++) {
			Map<SingleInterval, Double> vals=getCellProfile(chrX, numRNACopies);
			BEDFileIO.writeBEDGraph(vals, save+".random"+i+".bedgraph");
			Collection<Gene> genesBound=getBound(vals, genes);
			double ratio=(double)genesBound.size()/(double)genes.size();
			System.out.println(i+"\t"+genesBound.size()+"\t"+genes.size()+"\t"+ratio);
			
			String col="Cell"+i;
			for(Gene g: genesBound) {
				String name=g.toUCSC();
				genesByCell.set(name, col,1.0);
			}
			
		}
		
		genesByCell.write(save+".genesByCell.matrix");
		
	}

	private List<String> getColumns(int numPerms) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<numPerms; i++) {
			rtrn.add("Cell"+i);
		}
		
		return rtrn;
	}

	private List<String> getRow(IntervalTree<Gene> genes) {
		List<String> rtrn=new ArrayList<String>();
		
		Iterator<Gene> iter=genes.valueIterator();
		while(iter.hasNext()) {
			rtrn.add(iter.next().toUCSC());
		}
		
		
		return rtrn;
	}

	private Collection<Gene> getBound(Map<SingleInterval, Double> vals, IntervalTree<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		for(SingleInterval r: vals.keySet()) {
			Iterator<Gene> iter=genes.overlappingValueIterator(r.getReferenceStartPosition(), r.getReferenceEndPosition());
			while(iter.hasNext()) {
				rtrn.add(iter.next());
			}
		}
		
		return rtrn;
	}

	private Map<SingleInterval, Double> getCellProfile(SingleInterval region, int numRNACopies) {
		Map<SingleInterval, Double> rtrn=new TreeMap<SingleInterval, Double>();
		
		for(int i=0; i<numRNACopies; i++) {
			SingleInterval pos=samplePosition(region);
			double count=0;
			if (rtrn.containsKey(pos)) {count=rtrn.get(pos);}
			count++;
			rtrn.put(pos, count);
		}
		
		return rtrn;
	}

	private SingleInterval samplePosition(SingleInterval region) {
		int start=(int)Math.round(Math.random()*region.getReferenceEndPosition());
		int end=start+17000;
		return new SingleInterval(region.getReferenceName(), start, end);
	}
	
	public static void main(String[] args) throws IOException {
		String save=args[0];
		int numPerm=Integer.parseInt(args[1]);
		IntervalTree<Gene> genes=BEDFileIO.loadTree(args[2]).get("chrX");
		new SimulateXistLocalization(save, numPerm, genes);
	}
	
}
