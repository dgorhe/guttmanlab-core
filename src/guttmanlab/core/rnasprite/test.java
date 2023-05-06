package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class test {

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		RNAFilter filter=new RNAFilter();
		filter.addRNA("Xist");
		data.addFilter(filter);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			System.out.println(c.toString());
			counter++;
			if(counter%100==0){System.err.println(counter);}
		}
		data.close();
		System.err.println(counter);
	}

	private static List<String> getNames(Map<SingleInterval, Double> frequency) {
		List<String> rtrn=new ArrayList<String>();
		for(SingleInterval region: frequency.keySet()){
			rtrn.add(region.toUCSC());
		}
		return rtrn;
	}
	
}
