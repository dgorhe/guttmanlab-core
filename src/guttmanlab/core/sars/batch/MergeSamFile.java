package guttmanlab.core.sars.batch;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.datastructures.Pair;

public class MergeSamFile {
	
	
	private static void mergeSAM(File[] files, String output) throws IOException, InterruptedException{
		String cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar MergeSamFiles";
		
		for(int i=0; i<files.length; i++){
			cmd+=" I="+files[i].getName();
		}
		
		cmd+=" O="+output;
		cmd+=" SO=coordinate";
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name=merge --output stdout --error sterr"+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		String output=args[1];
		mergeSAM(files, output);
		
		
	}

	
	
}
