package guttmanlab.core.sars.batch;

import java.io.File;
import java.io.IOException;

public class BatchSubmitSecondRead {

	private static void secondRead(String input, String output, String name) throws IOException, InterruptedException{
		String cmd="java -jar /groups/guttman/mguttman/scripts/sars/PlotSecondReadPileups.jar "+input+" "+output;
				
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
		File[] files=new File(args[0]).listFiles();
		String saveDir=args[1];
		
		for(int i=0; i<files.length; i++){
			if(files[i].getAbsolutePath().endsWith(".sam") || files[i].getAbsolutePath().endsWith(".bam") ){
				String name=files[i].getName().replace(".sam", "");
				name=name.replace(".bam", "");
				String save=saveDir+"/"+name+".secondreads.bedgraph";
				secondRead(files[i].getAbsolutePath(), save, name);
			}
		}
		}
		else{System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=directory \n args[1]=save directory";
	
}
