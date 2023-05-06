package guttmanlab.core.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;

public class SubsetBedgraph {

	public SubsetBedgraph(Collection<SingleInterval> regions, Map<SingleInterval, Double> scores, String save) throws IOException {
		
		Map<String, IntervalTree<SingleInterval>> tree=ChIPUtils.makeTree(scores.keySet());
		
		Collection<SingleInterval> subSetList=new TreeSet<SingleInterval>();
		for(SingleInterval region: regions) {
			Collection<SingleInterval> overlappers=ChIPUtils.getRegions(tree, region);
			subSetList.addAll(overlappers);
		}
		
		FileWriter writer=new FileWriter(save);
		for(SingleInterval r: subSetList) {
			double score=scores.get(r);
			writer.write(r.toBedgraph(score)+"\n");
		}
		writer.close();
	}
	
	
	public static void main(String[] args) throws IOException {
		Collection<SingleInterval> regions=BEDFileIO.loadSingleIntervalFromFile(args[0]);
		Map<SingleInterval, Double> scores=BEDFileIO.loadbedgraph(new File(args[1]));
		String save=args[2];
		new SubsetBedgraph(regions, scores, save);
	}
	
	
}
