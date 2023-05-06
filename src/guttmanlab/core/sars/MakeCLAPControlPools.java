package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;

public class MakeCLAPControlPools {

	
	private static void mergeSAM(File[] files, String exclude, String name) throws IOException {
		String cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar MergeSamFiles";
		
		String save="control."+name+".bam";
		
		for(int i=0; i<files.length; i++){
			if(!files[i].getAbsolutePath().equals(exclude)){
				cmd+=" I="+files[i].getAbsolutePath();
			}
		}
		
		cmd+=" O="+save+" SO=coordinate";
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
		
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		
		for(int i=0; i<files.length; i++){
			String name=files[i].getName().replace(".bam", "");
			String exclude=files[i].getAbsolutePath();
			
			mergeSAM(files, exclude, name);
		}
		
		
	}

	
	
}
