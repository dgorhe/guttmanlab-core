package guttmanlab.core.util;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;

public class Join4 {
	
	private static void merge(String file1, String file2, String file3, String file4) throws IOException {
		Map<String, String> map1=parse(file1);
		Map<String, String> map2=parse(file2);
		Map<String, String> map3=parse(file3);
		Map<String, String> map4=parse(file4);
		
		Collection<String> allKeys=new TreeSet<String>();
		allKeys.addAll(map1.keySet());
		allKeys.addAll(map2.keySet());
		allKeys.addAll(map3.keySet());
		allKeys.addAll(map4.keySet());
		
		for(String key: allKeys) {
			String str1=get(map1, key);
			String str2=get(map2, key);
			String str3=get(map3, key);
			String str4=get(map4, key);
			System.out.println(str1+"\t"+str2+"\t"+str3+"\t"+str4);
		}
		
	}
	

	private static Map<String, String> parse(String file1) throws IOException {
		List<String> lines=BEDFileIO.loadLines(file1);
		
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String line: lines) {
			String key=line.split("\t")[0];
			rtrn.put(key, line);
		}
		
		return rtrn;
	}


	private static String get(Map<String, String> map1, String key) {
		if(map1.containsKey(key)) {return map1.get(key);}
		return key+"\tNA\tNA\tNA\tNA\tNA\tNA";
	}


	public static void main(String[] args) throws IOException {
		String file1= args[0];
		String file2= args[1];
		String file3= args[2];
		String file4= args[3];
		
		merge(file1, file2, file3, file4);
	}

	
	
}
