package guttmanlab.core.sars.batch;

import java.io.File;
import java.io.IOException;

public class BatchBin {

	private static void batchBin(String input, String output, String name, int binSize) throws IOException, InterruptedException{
		String cmd="java -jar /groups/guttman/mguttman/scripts/sars/BinCounts.jar "+input+" "+output+" "+binSize;
				
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
		File[] files=new File(args[0]).listFiles();
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith(".bam") || files[i].getName().endsWith(".sam")){
				String name=files[i].getName().replace(".bam", "");
				name=name.replace(".sam", "");
				String input=files[i].getAbsolutePath();
				int binSize=new Integer(args[2]);
				String output=args[1]+"/"+name+".w"+binSize;
				
				batchBin(input, output, name, binSize);
			}
		}
		
		}else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam files \n args[1]=save dir \n args[2]=bin size";
}
