package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.rnasprite.BarcodingDataStreaming;
import guttmanlab.core.rnasprite.HeatmapForCompartment;
import guttmanlab.core.rnasprite.SplitSymetricMatrix;
import guttmanlab.core.simulation.CoordinateSpace;

public class DNAContactMatrix {
	
	private static MatrixWithHeaders writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /groups/guttman/SPRITE2/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return new MatrixWithHeaders(rtrn);
	}
	
	private static void run(BarcodingDataStreaming data1, BarcodingDataStreaming data2, int resolution, String save, boolean weight) throws IOException, InterruptedException {
		MatrixWithHeaders mwh1=data1.getDNADNAContactMatrix(resolution, weight, CoordinateSpace.MM10.getRefSizes());
		MatrixWithHeaders mwh2=data2.getDNADNAContactMatrix(resolution, weight, CoordinateSpace.MM10.getRefSizes());
		
		
		MatrixWithHeaders norm1= writeAndNorm(mwh1, save+".d1");
		MatrixWithHeaders norm2= writeAndNorm(mwh2, save+".d2");
		
		
		norm1=HeatmapForCompartment.norm(norm1, resolution);
		norm2=HeatmapForCompartment.norm(norm2, resolution);
		
		
		diff(norm1, norm2, save+".diff");
		
		new SplitSymetricMatrix(norm1, norm2, save);
		//mwh.write(save);
		data1.close();
		data2.close();
	}
	
	private static void diff(MatrixWithHeaders norm1, MatrixWithHeaders norm2, String string) throws IOException {
		MatrixWithHeaders diff=new MatrixWithHeaders(norm1.getRowNames(), norm1.getColumnNames());
		for(String row: norm1.getRowNames()) {
			for(String column: norm1.getColumnNames()) {
				double val1=norm1.get(row, column);
				double val2=norm2.get(row, column);
				double diffVal=val2-val1;
				diff.set(row, column, diffVal);
			}
		}
		diff.write(string);
	}

	private static void run(BarcodingDataStreaming data1, BarcodingDataStreaming data2, int resolution, String save, boolean weight, SingleInterval region) throws IOException, InterruptedException {
		MatrixWithHeaders mwh1=data1.getDNADNAContactMatrix(region.getReferenceName(), resolution, weight, CoordinateSpace.MM10.getRefSizes().get(region.getReferenceName()));
		MatrixWithHeaders mwh2=data2.getDNADNAContactMatrix(region.getReferenceName(), resolution, weight, CoordinateSpace.MM10.getRefSizes().get(region.getReferenceName()));
		
		
		MatrixWithHeaders norm1= writeAndNorm(mwh1, save+".d1");
		MatrixWithHeaders norm2= writeAndNorm(mwh2, save+".d2");
		
		
		norm1=HeatmapForCompartment.norm(norm1, resolution, region);
		norm2=HeatmapForCompartment.norm(norm2, resolution, region);
		
		new SplitSymetricMatrix(norm1, norm2, save);
		//mwh.write(save);
		data1.close();
		data2.close();
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>4){
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			BarcodingDataStreaming data2=new BarcodingDataStreaming(new File(args[1]));
			int resolution=Integer.parseInt(args[2]);
			String save=args[3];
			boolean weight=Boolean.parseBoolean(args[4]);
			
			if(args.length>5) {
				SingleInterval region=new SingleInterval(args[5]);
				run(data, data2, resolution, save, weight, region);
			}
			else {
				run(data, data2, resolution, save, weight);
			}
			
			
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=clusters2 \n args[2]=resolution \n args[3]=save \n args[4]=weight \n args[5]=region (optional)";

}
