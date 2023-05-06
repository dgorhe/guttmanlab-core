package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class RenameGenes {

	public RenameGenes(BarcodingDataStreaming data, Map<String, String> alias, String rnaToPull, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Collection<String> allNames=getNames(rnaToPull, alias);
		
		int counter=0;
		while(data.hasNext()){
			Cluster c=data.next();
			boolean isGood=isCluster(c, allNames);
			if(isGood){
				Cluster renamed=c.renameRNA(alias);
				writer.write(renamed+"\n");
			}
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		
		data.close();
		
		writer.close();
	}
	
	private Collection<String> getNames(String rnaToPull, Map<String, String> alias) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String uid: alias.keySet()){
			String name=alias.get(uid);
			if(name.equalsIgnoreCase(rnaToPull)){rtrn.add(uid);}
		}
		
		return rtrn;
	}
	
	private static Map<String, String> parse(String string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			rtrn.put(line.split("\t")[0], line.split("\t")[1]);
		}
		
		return rtrn;
	}
	

	private boolean isCluster(Cluster c, Collection<String> allNames) {
		for(String name: allNames){
			if(c.getRNANames().contains(name)){return true;}
		}
		return false;
	}

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		Map<String, String> alias=parse(args[1]);
		String rna=args[2];
		String save=args[3];
		new RenameGenes(data, alias, rna, save);
	}

	
}
