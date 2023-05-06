package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class ContactMatrixAcrossChromosomes {
	
	private static File writeAndNorm(MatrixWithHeaders observed, String input) throws IOException, InterruptedException {
		observed.write(input);
		
		String cmd="python /Users/mguttman/Desktop/SPRITE_ActD/ice_matrix.py "+input;
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		File rtrn=new File(input+".iced");
		return rtrn;
	}
	
	
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String save=args[1];
		boolean weight=Boolean.parseBoolean(args[2]);
		
		List<String> chromosomes=new ArrayList<String>();
		chromosomes.add("chr1");
		chromosomes.add("chr2");
		chromosomes.add("chr3");
		chromosomes.add("chr4");
		chromosomes.add("chr5");
		chromosomes.add("chr6");
		chromosomes.add("chr7");
		chromosomes.add("chr8");
		chromosomes.add("chr9");
		chromosomes.add("chr10");
		chromosomes.add("chr11");
		chromosomes.add("chr12");
		chromosomes.add("chr13");
		chromosomes.add("chr14");
		chromosomes.add("chr15");
		chromosomes.add("chr16");
		chromosomes.add("chr17");
		chromosomes.add("chr18");
		chromosomes.add("chr19");
		chromosomes.add("chrX");
		
		MatrixWithHeaders mwh=data.getDNADNAContactMatrixByChromosome(weight, chromosomes);
		
		
		writeAndNorm(mwh, save);
		//mwh.write(save);
		data.close();
		
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=clusters \n args[1]=save \n args[2]=weight";

}
