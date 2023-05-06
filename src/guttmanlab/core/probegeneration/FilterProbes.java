package guttmanlab.core.probegeneration;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.annotation.io.BEDFileIO;

public class FilterProbes {

	private static boolean rejectProbe(String probe) {
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 15, 12);
		//RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe);
	}
	
	public static void main(String[] args) throws IOException{
		String probe="AAAAACTGCTGCTG";
		
		/*Collection<String> lines=BEDFileIO.loadLines(args[0]);
		FileWriter writer=new FileWriter(args[1]);
		int counter=0;
		for(String line: lines){
			String probe=line.split("\t")[7];
			boolean reject=rejectProbe(probe);
			if(!reject){writer.write(probe+"\n");}
			else{System.err.println(counter+" "+probe);}
			counter++;
		}
		writer.close();*/
	}
	
}
