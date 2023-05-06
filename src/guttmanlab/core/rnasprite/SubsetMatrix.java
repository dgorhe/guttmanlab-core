package guttmanlab.core.rnasprite;

import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

public class SubsetMatrix {

	
	public static void main(String[] args) throws IOException{
		List<String> lines=BEDFileIO.loadLines(args[0]);
		String geneName=args[1];
		String save=args[2];
		
		FileWriter writer=new FileWriter(save);
		
		String geneLine="";
		for(String line: lines){
			if(line.split("\t")[0].equals(geneName)){
				geneLine=line;
			}
		}
		
		String[] tokens=geneLine.split("\t");
		String[] headerTokens=lines.get(0).split("\t");
		for(int i=1; i<headerTokens.length; i++){
			if(!headerTokens[i].contains("Input")){
				System.err.println(headerTokens[i]);
				SingleInterval region=new SingleInterval(headerTokens[i]);
				writer.write(region.getReferenceName()+"\t"+region.getReferenceStartPosition()+"\t"+region.getReferenceEndPosition()+"\t"+tokens[i]+"\n");
			}
		}
		
		writer.close();
	}
	
}
