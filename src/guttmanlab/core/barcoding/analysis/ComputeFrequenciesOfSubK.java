package guttmanlab.core.barcoding.analysis;

import java.util.Collection;

public class ComputeFrequenciesOfSubK {

	public ComputeFrequenciesOfSubK(BarcodingData data, Cluster cluster){
		int score=data.quantify(cluster);
		
		Collection<Cluster> subK=cluster.enumerateKmers(cluster.getClusterSize()-1);
		
		//For ABC: Compute freq of C for all clusters containing AB relative to C versus all clusters
	}
	
}
