package guttmanlab.core.sars.batch;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.datastructures.Pair;

public class BatchMarkDuplicates {
	
	
	private static void sortSam(String input, String output, String name) throws IOException, InterruptedException{
		String cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar MarkDuplicates REMOVE_DUPLICATES=true M="+name+".metrics I="+input+" O="+output;
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith(".bam")){
				String name=files[i].getName().replace(".bam", "");
				String input=files[i].getAbsolutePath();
				String output=input.replace(".bam", ".noduplicates.bam");
				sortSam(input, output, name);
			}
		}
		
		
	}

	
	
}
