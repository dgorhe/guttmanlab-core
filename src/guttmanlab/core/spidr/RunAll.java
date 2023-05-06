package guttmanlab.core.spidr;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class RunAll {

	
	private static Map<File, Double> parseFiles(String input) throws IOException {
		Map<File, Double> rtrn=new TreeMap<File, Double>();
		
		List<String> lines=BEDFileIO.loadLines(input);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			rtrn.put(new File(tokens[0]), Double.parseDouble(tokens[1]));
		}
		
		return rtrn;
	}
	
	private static double sum(Map<File, Double> allFiles, Collection<File> controlFiles) {
		double sum=0;
		for(File f: controlFiles) {
			sum+=allFiles.get(f);
		}
		return sum;
	}
	
	
	private static Collection<File> getControls(File file, Map<File, String> groups) {
		Collection<File> rtrn=new TreeSet<File>();
		String group=groups.get(file);
		for(File f: groups.keySet()) {
			if(!groups.get(f).equals(group)) {rtrn.add(f);}
		}
		return rtrn;
	}
	
	private static void run(File file, Collection<File> controlFiles, double p, String saveDir, int numPerm, String binSizes) throws IOException {
		String name=file.getName().split("\\.")[0];
		String save=saveDir+"/"+name+".diff.bed";

		String cmd="java -jar -Xmx16000m /groups/guttman/mguttman/scripts/RandomDownsampleSignificant.jar";
		cmd+=" "+file.getAbsolutePath();
		cmd+=" "+format(controlFiles);
		cmd+=" "+p;
		cmd+=" "+numPerm;
		cmd+=" "+save;
		cmd+=" "+binSizes;
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);	
		System.out.println(cmd);
		
	}
	
	private static void runBigwig(File file, Collection<File> controlFiles, double p, String saveDir, int numPerm, String sizes) throws IOException {
		String name=file.getName().split("\\.")[0]+"_BW";
		String save=saveDir+"/"+name+".diff";

		String cmd="java -jar -Xmx16000m /groups/guttman/mguttman/scripts/RandomDownsample.jar";
		cmd+=" "+file.getAbsolutePath();
		cmd+=" "+format(controlFiles);
		cmd+=" "+p;
		cmd+=" "+numPerm;
		cmd+=" "+save;
		cmd+=" "+sizes;
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);	
		System.out.println(cmd);
		
	}
	
	private static String format(Collection<File> controlFiles) {
		String rtrn="";
		
		for(File f: controlFiles) {
			rtrn+=f.getAbsolutePath()+",";
		}
		
		return rtrn.substring(0, rtrn.length()-1);
	}
	
	private static Map<File, String> excludeFiles(String input) throws IOException {
		Map<File, String> rtrn=new TreeMap<File, String>();
		
		List<String> lines=BEDFileIO.loadLines(input);
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			File self=new File(tokens[0]);
			
			rtrn.put(self, tokens[2]);
		}
		
		return rtrn;
		
	}
	
	private static File[] subset(File[] allBams, int index) {
		File[] rtrn=new File[allBams.length-1];
		
		int counter=0;
		for(int i=0; i<allBams.length; i++) {
			if(i!=index) {
				rtrn[counter]=allBams[i];
				counter++;
			}
		}
		
		return rtrn;
	}

	
	private static void run(File sampleBam, File[] controlBams, String save) throws IOException {
		String name=sampleBam.getName().split("\\.")[0];
		
		String cmd="java -jar -Xmx16000m /groups/guttman/mguttman/scripts/RandomDownsample_CLIP.jar";
		cmd+=" "+sampleBam.getAbsolutePath();
		
		for(int i=0; i<controlBams.length; i++) {
			cmd+=" "+controlBams[i].getAbsolutePath();
		}
		cmd+=" "+save;
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);	
		System.out.println(cmd);
		
	}
	
	public static void main(String[] args) throws IOException {
		if(args.length>1) {
			
			File[] allBams=new File(args[0]).listFiles();
			String saveDir=args[1];
			
			for(int i=0; i<allBams.length; i++) {
				File sampleBam=allBams[i];
				File[] controlBAMs=subset(allBams, i);
				String save=saveDir+"/"+sampleBam.getName().split("\\.")[0];
				run(sampleBam, controlBAMs, save);
			}
		}
		else {System.err.println(usage);}
	}

	static String usage=" args[0]= all bams \n args[1]=save directory";

	

	

	
	
	
}
