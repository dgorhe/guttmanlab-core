package guttmanlab.core.splicing.speckle;

import java.io.File;
import java.io.IOException;

public class RunDedup {

	public static void main(String[] args) throws IOException, InterruptedException {
		File[] files=new File(args[0]).listFiles();
		String out=args[1];
		
		for(int i=0; i<files.length; i++) {
			String in=files[i].getAbsolutePath();
			String save=out+"/"+files[i].getName();
			String cmd="java -jar picard.jar MarkDuplicates -I "+in +" -M /Users/mguttman/Desktop/metrics -O "+save;
			System.err.println(cmd);
			Process p=Runtime.getRuntime().exec(cmd);
			p.waitFor();
		}
	}
	
}
