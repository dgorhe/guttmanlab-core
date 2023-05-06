package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.rnasprite.BarcodingDataStreaming;

public class SPRITEOverloadSimulation {

	public SPRITEOverloadSimulation(BarcodingDataStreaming data, int overloadRatio, String save) {
		
		
		
		
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			int overload=Integer.parseInt(args[1]);
			String save=args[2];
			new SPRITEOverloadSimulation(data, overload, save);
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=overload factor \n args[2]=save";
	
}
