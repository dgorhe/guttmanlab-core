package guttmanlab.core.clap.old;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class SigToBED {

	private static void write(Collection<SingleInterval> regions, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(SingleInterval region: regions){
			writer.write(region.toBED()+"\n");
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		if(args.length>1){
			double windowP=new Double(args[2]);
			double localP=new Double(args[3]);
			double enrichment=50.0;
			Collection<SingleInterval> sigWindows=BEDFileIO.loadCustom(args[0], windowP, localP, enrichment);
			String save=args[1];
			
			write(sigWindows, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=.sig format \n args[1]=save (.bed format) \n args[2]=window p \n args[3]=local p";
}
