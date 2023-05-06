package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class RankNormBedgraphs {

	private static void rankNorm(File[] inputFiles, String saveDir) throws IOException {
		Map<SingleInterval, Double> referenceValues=BEDFileIO.loadbedgraph(inputFiles[0]);
	}
	
}
