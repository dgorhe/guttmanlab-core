package guttmanlab.core.sharp;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;

public class AnnotatePeaks {

	public static void main(String[] args) throws IOException {
		if(args.length>3) {
			Map<SingleInterval, Double> peaks=BEDFileIO.loadbedgraph(new File(args[0]));
			Map<String, IntervalTree<String>> tree=BEDFileIO.loadGeneNamesFromRefFlat(args[1]);
			String save=args[2];
			double enrichment=new Double(args[3]).doubleValue();
			
			FileWriter writer=new FileWriter(save);
			for(SingleInterval peak: peaks.keySet()) {
				double score=peaks.get(peak);
				if(score>=enrichment) {
					Collection<String> names=find(peak, tree);
					writer.write(peak.toUCSC()+"\t"+score);
					for(String name: names) {writer.write("\t"+name);}
					writer.write("\n");
				}
			}
			writer.close();
		}
		else {System.err.println(usage);}
	}

	static String usage=" args[0]=peaks \n args[1]=ref flat \n args[2]=save \n args[3]=enrichment";
	
	private static Collection<String> find(SingleInterval peak, Map<String, IntervalTree<String>> tree) {
		Collection<String> rtrn=new TreeSet<String>();
		if(tree.containsKey(peak.getReferenceName())) {
			Iterator<String> iter=tree.get(peak.getReferenceName()).overlappingValueIterator(peak.getReferenceStartPosition(), peak.getReferenceEndPosition());
			while(iter.hasNext()) {
				String line=iter.next();
				rtrn.add(line.split("\t")[0]);
			}
		}
		return rtrn;
	}
	
	
}
