package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import guttmanlab.core.rnasprite.Cluster;

public class MergeIntoEnsemble {

	
	public MergeIntoEnsemble(File[] files, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(int i=0; i<files.length; i++) {
			System.err.println(files[i]);
			BarcodingDataStreaming data=new BarcodingDataStreaming(files[i]);
			while(data.hasNext()) {
				Cluster c=data.next();
				if(c.getClusterSize()>1) {
					writer.write(c.toDNAString()+"\n");
				}
			}
			data.close();
		}
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		new MergeIntoEnsemble(files, save);
	}
}
