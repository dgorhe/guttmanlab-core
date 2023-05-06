package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.DistanceToNuclearBody;
import guttmanlab.core.rnasprite.HeatmapForCompartment;
import guttmanlab.core.rnasprite.Kmer;
import guttmanlab.core.simulation.CoordinateSpace;

public class NormalizedDistanceToBody {
	
	private static MatrixWithHeaders writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return new MatrixWithHeaders(rtrn);
	}
	
	
	
	
	private static MatrixWithHeaders run(BarcodingDataStreaming data1, int resolution, String save, boolean weight) throws IOException, InterruptedException {
		MatrixWithHeaders mwh1=data1.getDNADNAContactMatrix(resolution, weight, CoordinateSpace.MM10.getRefSizes());
		//MatrixWithHeaders norm1= writeAndNorm(mwh1, save+".d1");
		MatrixWithHeaders norm1=HeatmapForCompartment.norm(mwh1, resolution);
		norm1.write(save+".matrix");
		data1.close();
		return norm1;
		
		
		//new SplitSymetricMatrix(norm1, mwh1, save+".split.matrix");
	}
	
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			int resolution=Integer.parseInt(args[1]);
			String save=args[2];
			boolean weight=Boolean.parseBoolean(args[3]);
			
			Kmer kmer=new Kmer();
			kmer.addIntervals(BEDFileIO.loadSingleIntervalFromFile(args[4]));
			
			MatrixWithHeaders norm1=run(data, resolution, save, weight);
			//new DistanceToNuclearBody(norm1, kmer, save+".bedgraph");
			
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=resolution \n args[2]=save \n args[3]=weight \n args[4]=compartment regions";

}
