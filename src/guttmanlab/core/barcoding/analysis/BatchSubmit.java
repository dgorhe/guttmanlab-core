package guttmanlab.core.barcoding.analysis;

import java.io.File;
import java.io.IOException;

public class BatchSubmit {

	public static void main(String[] args) throws IOException, InterruptedException{
		File[] files=new File(args[0]).listFiles();
		String saveDir=args[1];
		int binSize=new Integer(args[2]);
		int k=new Integer(args[3]);
		int maxCluster=new Integer(args[4]);
		
		for(int i=0; i<files.length; i++){
			String jobName="Enumerate"+i;
			
			File file=files[i];
			String save=saveDir+"/"+file.getName();
			String command="java8 -jar /storage/mguttman/Software/EnumerateKMers.jar "+file.getAbsolutePath()+" "+save+" "+binSize+" "+k+" "+maxCluster;
			String exec="qsub -b y -N "+jobName+" -cwd -V -m bae -M mguttman@caltech.edu ";
			Process p=Runtime.getRuntime().exec(exec+command);
			p.waitFor();
			System.err.println(exec+command);
		}
		
	}
	
}
