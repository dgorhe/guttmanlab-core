package guttmanlab.core.chip;

import java.io.File;
import java.io.IOException;
import java.util.List;

public class BatchChIP {

	
	/*private static void bin(File file, File control, String outdir) throws IOException {
		String cmd="java -jar -Xmx6000m /groups/guttman/mguttman/scripts/ChIPUtils.jar ";
		
		
		String save=outdir+"/"+file.getName()+".enrich";
		cmd+=file.getAbsolutePath()+" "+control.getAbsolutePath()+" "+save;
		
		String name=file.getName();
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output "+outdir+"/"+name+".out"+" --error "+outdir+"/"+name+".err"+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		
		System.err.println(cmd);
		
	}*/
	
	
	private static void run(File bam, String save, int windowSize) throws IOException {
		String cmd="java -jar /groups/guttman/mguttman/scripts/ChIPUtils.jar "+bam.getAbsolutePath()+" "+save+" "+windowSize; 
		String name=bam.getName().split("\\.")[0];
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		System.err.println(cmd);
	}
	
	
	private static void mergeAllExcept(File[] files, int pos, String saveDir) throws IOException {
		String name=files[pos].getName().split("\\.")[0];
		String cmd="java -jar /groups/guttman/software/picard.2.18.7/picard.jar MergeSamFiles O="+saveDir+"/"+name+".control.bam SO=coordinate";
		
		for(int i=0; i<files.length; i++) {
			if(files[i].getName().endsWith(".bam") && i!=pos) {
				cmd=cmd+" I="+files[i].getAbsolutePath();
			}
		}
		
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		System.err.println(cmd);
		
	}

	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>1) {
		File[] files=new File(args[0]).listFiles();
		System.err.println(args[0]+" "+files.length);
		
		String outdir=args[1];
		int windowSize=Integer.parseInt(args[2]);
		
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith(".bam")){
				String save=outdir+"/"+files[i].getName()+".w"+windowSize+".bed";
				run(files[i], save, windowSize);;
				//run(files[i], save, windowSize);
			}
			
			
		}
		
		}else {System.err.println(usage);}
	}

	
	

	static String usage=" args[0]=input dir \n args[1]=save dir";
	
	
}
