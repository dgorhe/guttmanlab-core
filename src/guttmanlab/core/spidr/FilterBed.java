package guttmanlab.core.spidr;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class FilterBed {

	public static void main(String[] args) throws IOException {
		Map<SingleInterval, Double> scores=BEDFileIO.loadbedgraph(new File(args[0]));
		double cutoff=Double.parseDouble(args[1]);
		
		for(SingleInterval region: scores.keySet()) {
			double score=scores.get(region);
			if(score>cutoff) {System.out.println(region.toBedgraph(score));}
		}
	}
	
}
