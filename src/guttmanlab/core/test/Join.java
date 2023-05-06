package guttmanlab.core.test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class Join {

	public Join(File file1, File file2, String save) throws IOException{
		Map<String, String> lines1=parse(file1);
		Map<String, String> lines2=parse(file2);
		
		FileWriter writer=new FileWriter(save);
		
		for(String key: lines1.keySet()){
			if(lines2.containsKey(key)){
				writer.write(key+"\t"+lines1.get(key)+"\t"+lines2.get(key)+"\n");
			}
		}
		
		writer.close();
	}
	
	private Map<String, String> parse(File file1) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		Collection<String> lines=BEDFileIO.loadLines(file1.getAbsolutePath());
		
		for(String line: lines){
			String key=line.split("\t")[0];
			rtrn.put(key, line);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		File file1=new File(args[0]);
		File file2=new File(args[1]);
		String save=args[2];
		new Join(file1, file2, save);
	}
}
