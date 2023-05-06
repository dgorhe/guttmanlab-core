package guttmanlab.core.clap.old;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Set;
import java.util.TreeSet;

public class MakeUnique {

	public MakeUnique(File file, String save) throws IOException{
		Set<String> lines=readLines(file);
		write(save, lines);
	}

	private void write(String save, Set<String> lines) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String line: lines){
			writer.write(line+"\n");
		}
		
		writer.close();
	}

	private Set<String> readLines(File file) throws IOException {
		Set<String> set=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine=reader.readLine();
		
		while((nextLine=reader.readLine())!=null){
			set.add(nextLine);
		}
		
		reader.close();
		return set;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File file=new File(args[0]);
			String save=args[1];
			new MakeUnique(file, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage="args[0]=input \nargs[1]=output";
}
