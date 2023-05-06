package guttmanlab.core.splicing.speckle;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class GetGenesWithinClusters {

	public static Collection<String> getGenesWithinClusters(Collection<SingleInterval> regions, Map<String, IntervalTree<SingleInterval>> geneTree) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(SingleInterval region: regions) {
			Collection<SingleInterval> overlappingGenes=getOverlaps(region, geneTree);
			for(SingleInterval r: overlappingGenes) {
				rtrn.add(r.getName());
			}
		}
		return rtrn;
	}

	private static Collection<SingleInterval> getOverlaps(SingleInterval region, Map<String, IntervalTree<SingleInterval>> geneTree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		String chr=region.getReferenceName();
		
		if(geneTree.containsKey(chr)) {
			IntervalTree<SingleInterval> tree=geneTree.get(chr);
			Iterator<SingleInterval> iter=tree.overlappingValueIterator(region.getReferenceStartPosition(), region.getReferenceEndPosition());
			while(iter.hasNext()) {
				rtrn.add(iter.next());
			}
		}
		return rtrn;
	}
	
	private static MatrixWithHeaders makeMatrix(File[] files) throws IOException {
		List<String> columns=getNames(files);
		List<String> rows=getRows(files);
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i<files.length; i++) {
			String col=files[i].getName();
			List<String> vals=BEDFileIO.loadLines(files[i].getAbsolutePath());
			for(String row: vals) {
				rtrn.set(row, col, 1);
			}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalsFromFile(args[0]);
		Map<String, IntervalTree<SingleInterval>> geneTree=BEDFileIO.loadSingleIntervalTreeFromGTF(args[1], true);
		Collection<String> genes=getGenesWithinClusters(regions, geneTree);
		for(String g: genes) {System.out.println(g);}
	}
	
	private static List<String> getRows(File[] files) throws IOException {
		Collection<String> set=new TreeSet<String>();
		for(int i=0; i<files.length; i++) {
			set.addAll(BEDFileIO.loadLines(files[i].getAbsolutePath()));
		}
		List<String> rtrn=new ArrayList<String>();
		rtrn.addAll(set);
		return rtrn;
	}

	private static List<String> getNames(File[] files) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<files.length; i++) {
			rtrn.add(files[i].getName());
		}
		return rtrn;
	}

	/*public static void main(String[] args) throws IOException {
		File[] files=new File(args[0]).listFiles();
		MatrixWithHeaders matrix= makeMatrix(files);
		matrix.write(args[1]);
	}*/

	
	
}
