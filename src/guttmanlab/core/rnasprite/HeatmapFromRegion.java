package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class HeatmapFromRegion {

	
	public static void main(String[] args) throws IOException{
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			SingleInterval region=new SingleInterval(args[1]);
			int binResolution=new Integer(args[2]);
			boolean weight=new Boolean(args[3]);
			String save=args[4];
			MatrixWithHeaders mwh=data.getDNADNAContactMatrix(binResolution, region, weight);
			mwh.write(save);
			data.close();
			
			System.out.println("ID\tscore");
			for(String row: mwh.getRowNames()){
				double sum=mwh.get(row, row);
				System.out.println(row+"\t"+sum);
			}
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=region (chr:start-end) \n args[2]=bin resolution \n args[3]=weight by cluster size (true or false) \n args[4]=save";
}
