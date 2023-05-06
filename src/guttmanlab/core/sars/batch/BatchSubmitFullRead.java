package guttmanlab.core.sars.batch;

import java.io.File;
import java.io.IOException;

public class BatchSubmitFullRead {

	private static void secondRead(String input, String output, String name, int binSize) throws IOException, InterruptedException{
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
		String saveDir=args[1];
		int binSize=new Integer(args[2]);
		
		for(int i=0; i<files.length; i++){
			if(files[i].getAbsolutePath().endsWith(".sam") || files[i].getAbsolutePath().endsWith(".bam") ){
				String name=files[i].getName().replace(".sam", "");
				name=name.replace(".bam", "");
				String save=saveDir+"/"+name+".w"+binSize+"fullreads";
				secondRead(files[i].getAbsolutePath(), save, name, binSize);
			}
		}
		}
		else{System.err.println(usage);}
		
	}
	
	static String usage=" args[0]=directory \n args[1]=save directory \n args[2]=bin size";
	
}
