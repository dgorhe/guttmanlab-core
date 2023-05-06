package guttmanlab.core.sequence;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;


public class SplitFasta {

	public static void main(String[] args) throws IOException {
		if(args.length>4) {
		Collection<Sequence> sequences=FastaFileIOImpl.readFromFile(args[0]);
		String dir=args[1];
		int fileSize=Integer.parseInt(args[2]);
		
		FileWriter writer=new FileWriter(dir+"/0.fa");
		int counter=0;
		for(Sequence seq: sequences) {
			//String name=seq.getName().replace("_intron", "");
			//name=name.split(":")[1];
			//seq.setName(name);
			writer.write(seq.toFasta()+"\n");
			counter++;
			if(counter%fileSize ==0) {
				int fileNum=counter/fileSize;
				writer.close();
				writer=new FileWriter(dir+"/"+fileNum+".fa");
			}
		}
		writer.close();
		
		
		
		File[] files=new File(dir).listFiles();
		for(int i=0; i<files.length; i++) {
			run(files[i], args[3], args[4]);
			
		}
		
		}
		else {System.err.println(usage);}
	}
	
	static String usage="args[0]= fasta \n args[1]=out dir \n args[2]=file size \n args[3]=gene bed \n args[4]=type";
	

	private static void run(File file, String geneBED, String type) throws IOException {
		String out=file.getAbsolutePath()+".primers";
		String cmd="java -jar DesignRTPrimers2.jar "+file.getAbsoluteFile()+" "+geneBED+" "+out+" "+type; 
		String name=file.getName();
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
		
		Runtime.getRuntime().exec(cmd);
		System.err.println(cmd);
	}
	
}
