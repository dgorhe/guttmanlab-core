package guttmanlab.core.sars;


import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class GetNSP1 {

	int n=100;
	public GetNSP1(String fasta, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Collection<Sequence> seqs=FastaFileIOImpl.readFromFile(fasta);
		for(Sequence seq: seqs){
			String first=seq.getSequenceBases().substring(0, n);
			
			
			
			writer.write(">"+seq.getName()+"\n"+first+"\n");
		}
		
		writer.close();
	}
	
	


	public static void main(String[] args) throws IOException{
		String fasta=args[0];
		String save=args[1];
		new GetNSP1(fasta, save);
	}
	
}
