package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class RNARNAHeatmap {

	
	
	public static void main(String[] args) throws IOException {
		if(args.length>2) {
			BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
			List<String> RNAList=parse(args[1]);
			MatrixWithHeaders mwh= data.getRNARNAContactMatrix(RNAList, false);
			mwh.write(args[2]);
		}
		else {System.err.println(usage);}
		
	}
	
	private static List<String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string);
		List<String> rtrn=new ArrayList<String>();
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.add(tokens[0]);
		}
		return rtrn;
	}

	static String usage=" args[0]=data \n args[1]=RNA list \n args[2]=save";
}
