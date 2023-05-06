package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.IOException;

public class BatchCounts {

	private static void batchCount(String cluster, String region, String binResolution, String save, String name) throws IOException, InterruptedException{
		String cmd="java -jar -Xmx8000m /groups/guttman/mguttman/scripts/SPRITE/ContactsByChr.jar "+cluster+" "+region+" "+binResolution+" "+save;
				
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		String region=args[1];
		String binResolution=args[2];
		String saveDir=args[3];
		
		for(int i=0; i<files.length; i++){
			String name=files[i].getName();
			String save=saveDir+"/"+name;
			batchCount(files[i].getAbsolutePath(), region, binResolution, save, name);
		}
	}
	
}
