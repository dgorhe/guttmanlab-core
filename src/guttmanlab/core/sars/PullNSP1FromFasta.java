package guttmanlab.core.sars;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class PullNSP1FromFasta {

	public PullNSP1FromFasta(File fasta, String save) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(fasta)));
		FileWriter writer=new FileWriter(save);
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			String name=nextLine;
			nextLine = reader.readLine();
			String seq=nextLine;
			if(name.startsWith(">NSP1|")){
				System.err.println(name);
				writer.write(name+"\n"+seq+"\n");
			}
			
		}
		reader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException{
		File fasta=new File(args[0]);
		String save=args[1];
		new PullNSP1FromFasta(fasta, save);
	}
	
}
