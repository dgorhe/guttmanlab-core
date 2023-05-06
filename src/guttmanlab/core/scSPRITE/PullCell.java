package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.Cluster;

public class PullCell {


	private static MatrixWithHeaders subMatrix(MatrixWithHeaders data, SingleInterval region) {
		
		List<String> list=new ArrayList<String>();
		for(String row: data.getRowNames()){
			SingleInterval r1=new SingleInterval(row);
			if(r1.overlaps(region)){list.add(row);}
		}
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(list, list);
		
		for(String r1: list){
			for(String r2: list){
				rtrn.set(r1, r2, data.get(r1, r2));
			}
		}
		return rtrn;
	}
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
			MatrixWithHeaders data=new MatrixWithHeaders(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			MatrixWithHeaders sub=subMatrix(data, region);
			sub.write(args[2]);
			
			/*String row=args[1];
			String column=args[2];
			double val=data.get(row, column);
			System.out.println(row+"\t"+column+"\t"+val);*/
		}
		else{System.err.println(usage);}
	}
	
	
	

	




	static String usage=" args[0]=matrix \n args[1]=region \n args[2]=save";
}
