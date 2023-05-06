package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class Translate {
	
	public Translate(File fasta, String save) throws IOException {
		
		FileWriter writer=new FileWriter(save);
		FastaFileIOImpl f=new FastaFileIOImpl(fasta);
	
		int counter=0;
		while(f.hasNext()) {
			Sequence seq=f.next();
			String protein=seq.translate();
			
			writer.write(">"+seq.getName()+"\n"+protein+"\n");
			counter++;
			if(counter%10000==0) {System.err.println(counter+" "+seq.getName()+" "+protein);}
		}
		writer.close();
		
	}
	
	public static void main(String[] args) throws IOException {
		File fasta=new File(args[0]);
		String save=args[1];
		new Translate(fasta, save);
	}

}
