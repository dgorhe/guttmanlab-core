package guttmanlab.core.barcodeidentification;

import java.io.File;
import java.io.IOException;

public class BatchSubmit {

	private static String jar="/groups/guttman/mguttman/scripts/PooledCLIP.jar";
	private static String allBarcodes="/groups/guttman/mguttman/scripts/barcodes.fa";
	private static String R8="/groups/guttman/mguttman/scripts/R8.fa";
	private static String R7="/groups/guttman/mguttman/scripts/R7.fa";
	private static String R6="/groups/guttman/mguttman/scripts/R6.fa";
	private static String R5="/groups/guttman/mguttman/scripts/R5.fa";
	private static String R4="/groups/guttman/mguttman/scripts/R4.fa";
	private static String R3="/groups/guttman/mguttman/scripts/R3.fa";
	private static String R2="/groups/guttman/mguttman/scripts/R2.fa";
	private static String pids="/groups/guttman/mguttman/scripts/pids.fa";
	
	private static void batchSubmit(String read1, String read2, String saveDir, String name) throws IOException, InterruptedException{
		String output=saveDir+"/"+name;
		
		String cmd="java -jar "+jar +" "+read1+" "+read2+" "+allBarcodes+" "+R8+" "+R7+" "+R6+" "+R5+" "+R4+" "+R3 +" "+R2+" "+pids+" "+output;
			
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	private static String getRead1(File dir) {
		for(File file: dir.listFiles()){
			if(file.getName().endsWith("R1.fastq.gz")){return file.getAbsolutePath();}
		}
		return null;
	}
	
	private static String getRead2(File dir) {
		for(File file: dir.listFiles()){
			if(file.getName().endsWith("R2.fastq.gz")){return file.getAbsolutePath();}
		}
		return null;
	}
	
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>1){
			File[] dirs=new File(args[0]).listFiles();
			String saveDir=args[1];
			
			for(File dir: dirs){
				if(dir.isDirectory()){
					String read1=getRead1(dir);
					String read2=getRead2(dir);
					String name=dir.getName();
					//System.err.println(dir.getAbsolutePath()+" "+name+" "+read1+" "+read2);
					
					if(read1!=null && read2!=null){
						batchSubmit(read1, read2, saveDir, name);
					}
				}
			}
			
		}
		else{System.err.println(usage);}
	}
	
	

	static String usage=" args[0]=input directory \n args[1]=save dir";
	
}
