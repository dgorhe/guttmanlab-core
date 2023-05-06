package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

public class PermuteRNAClusters {

	
	public PermuteRNAClusters(BarcodingDataStreaming data, String save) throws IOException {
		
		BarcodingDataStreaming shuffled=data.shuffleRNA(new File(save));
		
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>1) {
			
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			String save=args[1];
			new PermuteRNAClusters(data, save);
			
		}
		else {System.err.println(usage);}	
	}
	
	private static String usage=" args[0]=clusters \n args[1]=save";
	
}
