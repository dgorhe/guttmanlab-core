package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.simulation.CoordinateSpace;

public class DNAContactMatrix {
	
	private static File writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return rtrn;
	}
	
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>3){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		int resolution=Integer.parseInt(args[1]);
		String save=args[2];
		boolean weight=Boolean.parseBoolean(args[3]);
		String chr=null;
		
		if(args.length>4){
			chr=args[4];
		}
		
		MatrixWithHeaders mwh=null;
		if(chr!=null){
			mwh=data.getDNADNAContactMatrix(chr, resolution, weight, CoordinateSpace.MM10.getRefSizes().get(chr));
		}
		else{
			
			mwh=data.getDNADNAContactMatrix(resolution, weight, CoordinateSpace.MM10.getRefSizes());
		}
		
		
		writeAndNorm(mwh, save);
		//mwh.write(save);
		data.close();
		
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=resolution \n args[2]=save \n args[3]=weight \n (args[4]=optional region)";

}
