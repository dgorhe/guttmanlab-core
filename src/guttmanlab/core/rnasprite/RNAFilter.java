package guttmanlab.core.rnasprite;

import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.rnasprite.BarcodingDataStreaming.ClusterFilter;

public class RNAFilter implements ClusterFilter{
	Collection<String> rnas;
	
	public RNAFilter(Collection<String> rnas){
		this.rnas=rnas;
	}
	
	public RNAFilter(){
		this.rnas=new TreeSet<String>();
	}
	
	public void addRNA(String rna){this.rnas.add(rna);}
	
	@Override
	public boolean evaluate(Cluster c) {
		for(String rna: c.getRNANames()){
			if(rnas.contains(rna)){return true;}
		}
		return false;
	}

	
	
}
