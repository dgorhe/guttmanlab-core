package guttmanlab.core.clap.updated;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.sharp.AssignReads;

public class CLAPAnalysis_MainMethods {

	public static void assignReads(File bam, Map<String, IntervalTree<Annotation>> genes, String save) throws IOException {
		new AssignReads(bam, genes, save);
	}
	
	
}
