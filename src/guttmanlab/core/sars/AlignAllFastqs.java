package guttmanlab.core.sars;

import java.io.File;
import java.io.IOException;

import guttmanlab.core.datastructures.Pair;

public class AlignAllFastqs {

	private static void bowtie(String fq1, String fq2, String save, String name, String index, boolean local, boolean noUnalign) throws IOException, InterruptedException{
		
		String cmd="bowtie2";
		//System.err.println(cmd);
		
		if(local){
			cmd+=" --local";
		}
		
		if(noUnalign){
			cmd+=" --no-unal";
		}
		
		cmd+=" -x "+index+" -1 "+fq1+" -2 "+fq2 +" -S "+save;
		
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		//p.waitFor();
		//System.err.print(p.getErrorStream());
		
		System.err.println(cmd);
	}
	
	private static Pair<File> getReads(File[] files) {
		Pair<File> rtrn=new Pair<File>();
		
		for(int i=0; i<files.length; i++){
			if(!files[i].getName().startsWith("Trimmed")) {
				if(files[i].getName().endsWith("_R1.fastq.gz")){rtrn.setValue1(files[i]);}
				if(files[i].getName().endsWith("_R2.fastq.gz")) {rtrn.setValue2(files[i]);}
			}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException{
		if(args.length>2){
			//String path="/groups/guttman/primarydata/sequencingruns/Freeman_DM_CLAPpool_iLabs_14668/";
			String path=args[0];
			String index=args[1];
			String saveDir=args[2];
			boolean local=false;
			if(args.length>3){local=new Boolean(args[3]);}
			boolean noUnalign=false;
			if(args.length>4){noUnalign=new Boolean(args[4]);}
			
			File[] dir=new File(path).listFiles();
			for(int i=0; i<dir.length; i++){
				File currentDirectory=dir[i];
				if(currentDirectory.isDirectory()){
					File[] files=currentDirectory.listFiles();
					Pair<File> fq=getReads(files);
					String name=currentDirectory.getName();
					//System.err.println(name);
					String save=saveDir+"/"+name+".sam";
					if(fq.hasValue1() && fq.hasValue2()){
						bowtie(fq.getValue1().getAbsolutePath(), fq.getValue2().getAbsolutePath(), save, name, index, local, noUnalign);
					}
				}
			}
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=path \n args[1]=index \n args[2]=save directory \n args[3]=local alignment \n args[4]=no unaligned";
	
}
