package guttmanlab.core.snps;

import java.io.File;
import java.io.IOException;

public class BatchSplit {

	private static void batchBin(String bam, String snp, String save, String bed, String name) throws IOException, InterruptedException{
		String cmd="java -jar -Xmx8000m /groups/guttman/mguttman/scripts/SHARP/SplitBAMBySNP.jar "+bam+" "+snp+" "+bed+" "+save;
				
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	
	private static String getBam(File dir) {
		File[] files=dir.listFiles();
		
		for(int i=0; i<files.length; i++){
			if(files[i].getName().endsWith(".bam")){return files[i].getAbsolutePath();}
		}
		return null;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		String snp=args[0];
		String saveDir=args[1];
		String bed=args[2];		
		
		String dir="/groups/guttman/jamie/data/sharp/heard_data_nature_2020/correctFq/alignments/mm10_default/star/";
		File[] dirs=new File(dir).listFiles();
		for(int i=0; i<dirs.length; i++){
			String name=dirs[i].getName();
			String bam=getBam(dirs[i]);
			if(bam!=null){
				String save=saveDir+"/"+name+".split";
				batchBin(bam, snp, save, bed, name);
			}
		}
	}	
			
		
		
		
	
	
	

	static String usage=" args[0]=snp file \n args[1]=save dir \n args[2]=gene bed file";
}
