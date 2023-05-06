package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.datastructures.Pair;

public class SortSam {
	
	
	private static void sortSam(String input, String output, String name) throws IOException, InterruptedException{
		String cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar SortSam I="+input+" O="+output+" SO=coordinate";
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		
		for(int i=0; i<files.length; i++){
			String name=files[i].getName().replace(".sam", "");
			String input=files[i].getAbsolutePath();
			String output=input.replace(".sam", ".bam");
			sortSam(input, output, name);
		}
		
		
	}

	
	
}
