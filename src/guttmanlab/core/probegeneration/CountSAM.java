package guttmanlab.core.probegeneration;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class CountSAM {

	public static void main(String[] args) throws IOException {
		Map<String, Integer> counts=new TreeMap<String, Integer>();
		List<String> lines=BEDFileIO.loadLines(args[1]);
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(args[0])));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			
			String[] tokens=nextLine.split("\t");
			if(tokens.length>8) {
				String key=tokens[9];
				int count=0;
				if(counts.containsKey(key)) {
					count=counts.get(key);
				}
				count++;
				if(count>1) {System.err.println(key+" "+count);}
				counts.put(key, count);
			}
		}
		
		for(String line: lines) {
			String key=line.split("\t")[0];
			if(counts.containsKey(key)) {
				int count=counts.get(key);
				if(count==1) {System.out.println(line);}
				else {System.err.println(line+" "+count);}
			}
		}
		
		reader.close();
	}
	
}
