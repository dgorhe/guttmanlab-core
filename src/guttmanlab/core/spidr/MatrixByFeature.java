package guttmanlab.core.spidr;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.splicing.speckle.GTFToJunctions;

public class MatrixByFeature {

	public MatrixByFeature(GTFToJunctions genes) {
		Collection<SingleInterval> junctions=genes.getIntrons();
		
	}
	
	public static void main(String[] args) throws IOException {
		GTFToJunctions gtf=new GTFToJunctions(new File(args[0]));
		
		
		
		
		new MatrixByFeature(gtf);
	}
	
}
