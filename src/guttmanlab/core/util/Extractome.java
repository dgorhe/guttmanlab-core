package guttmanlab.core.util;

import java.io.File;
import java.io.IOException;

public class Extractome {

	private static String runExtractome(String regions, String genome, String saveDir, String name) throws IOException, InterruptedException {
		String cmd="/central/home/mguttman/yes/bin/python /groups/guttman/mguttman/scripts/IGVExtractome/extractome/extractome/extract.py ";
		
		cmd+=regions;
		cmd+=" --genome "+genome;
		cmd+=" --name "+name;
		cmd+=" --output "+saveDir;
		
		System.err.println(cmd);
		Process p=Runtime.getRuntime().exec(cmd);
		p.waitFor();
		
		return saveDir+"/"+name+".chain";
	}
	
	
	private static void runCrossMap(File[] files, String chain, String type, String saveDir) throws IOException, InterruptedException {
		
		
		String cmd="/central/home/mguttman/yes/bin/CrossMap.py " +type;
		cmd+=" "+chain;
		
		for(File file: files) {
			String save=saveDir+"/"+file.getName()+".converted"+getExtension(file);
			String name=file.getName().split("\\.")[0];
			String newCmd=cmd+" "+file.getAbsolutePath()+" "+save;
			//String qsub="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+newCmd+" &";
			System.err.println(newCmd);
			Process p=Runtime.getRuntime().exec(newCmd);
			p.waitFor();
		}
		
	}
	
	private static String getExtension(File file) {
		String[] tokens=file.getName().split("\\.");
		return "."+tokens[tokens.length-1];
	}


	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length>2) {
			String regions=args[0];
			String genome=args[1];
			String saveDir=args[2];
			String name=args[3];
			String chain=runExtractome(regions, genome, saveDir, name);
			/*String chain=args[0];
			File[] files=new File(args[1]).listFiles();
			String saveDir=args[2];
			String type="bed";
			runCrossMap(files, chain, type, saveDir);*/
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=regions to extract \n args[1]=genome (mm10, hg38) \n args[2]=save directory \n args[3]=name \n args[4]=files to convert";
	//static String usage=" args[0]=chain \n args[1]=files to convert \n args[2]=save directory";
	
}
