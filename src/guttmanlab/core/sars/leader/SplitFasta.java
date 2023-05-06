package guttmanlab.core.sars.leader;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class SplitFasta {

	public SplitFasta(String fasta, String save, String name) throws IOException {
		FileWriter writer=new FileWriter(save);
		Collection<Sequence> seqs=FastaFileIOImpl.readFromFile(fasta);
		
		for(Sequence seq: seqs) {
			if(seq.getName().startsWith(name)) {
				writer.write(seq.toFasta()+"\n");
			}
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException {
		String fasta=args[0];
		String save=args[1];
		String name=args[2];
		new SplitFasta(fasta, save, name);
	}
	
}
