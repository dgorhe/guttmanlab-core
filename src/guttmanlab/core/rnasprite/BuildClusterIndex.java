package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;

/**
 * Build a index to rapidly find all clusters containing a given RNA or genomic bin
 * @author mguttman
 *
 */
public class BuildClusterIndex {
	
	public BuildClusterIndex(File file, String save) throws IOException{
		
		FileWriter writer=new FileWriter(save);
		
		RandomAccessFile reader=new RandomAccessFile(file, "r");
		
		/*reader.seek(8158770);
		byte[] b=new byte[283];
		reader.read(b);
		String s=new String(b);
		System.err.println(s);
		
		reader.close();*/
		
		Map<String, List<long[]>> map=new TreeMap<String, List<long[]>>();
		
		long currentIndex=reader.getFilePointer();
		String nextLine=reader.readLine();
		while(nextLine!=null){
			Cluster c=parseCluster(nextLine);
			long[] array=new long[2];
			array[0]=currentIndex;
			array[1]=reader.getFilePointer()-currentIndex;
			
			for(String rna:c.getRNANames()){
				List<long[]> list=new ArrayList<long[]>();
				if(map.containsKey(rna)){list=map.get(rna);}
				list.add(array);
				map.put(rna, list);
			}
			
			
			
			
			
			currentIndex=reader.getFilePointer();
			
			nextLine=reader.readLine();
		}
		
		reader.close();
		
		
		for(String rna: map.keySet()){
			writer.write(rna);
			List<long[]> indexes=map.get(rna);
			for(long[] array: indexes){
				writer.write("\t"+array[0]+","+array[1]);
			}
			writer.write("\n");
		}
		
		writer.close();
		
		
	}

	private Cluster parseCluster(String nextLine) {
		String[] tokens=nextLine.split("\t");
		String barcode=tokens[0];
		
		Cluster cluster=new Cluster(barcode);
		cluster.setLine(nextLine);
		for(int i=1; i<tokens.length; i++){
			if(tokens[i].contains(":")){
				String name=tokens[i].split("_")[0];
				boolean dna=isDNA(name);
				tokens[i]=tokens[i].split("\\)_")[1];
				String chr=tokens[i].split(":")[0];
				int start=new Integer(tokens[i].split(":")[1].split("-")[0]);
				int end=new Integer(tokens[i].split(":")[1].split("-")[1]);
				
				
				
				
				SingleInterval interval=new SingleInterval(chr, start, end);
				if(dna){cluster.addDNARead(interval);}
				else{
					String rnaName=getName(name);
					RNAInterval i2=new RNAInterval(interval);
					i2.setName(rnaName);
					cluster.addRNARead(i2);
				}
			}
		}
		return cluster;
	}
	
	private boolean isDNA(String string) {
		if(string.startsWith("DPM")){return true;}
		return false;
	}
	
	private String getName(String string) {
		String name=string.split("\\(")[1].replaceAll("\\)", "");
		return name;
	}
	
	private static void readFromIndex(File file, File index, String gene) throws IOException{
		Map<String, List<long[]>> indexKeys=parse(index);
		List<long[]> list=indexKeys.get(gene);
		
		RandomAccessFile reader=new RandomAccessFile(file, "r");
		
		for(long[] array: list){
			reader.seek(array[0]);
			byte[] b=new byte[(int) array[1]-1];
			reader.read(b);
			String s=new String(b);
			System.err.println(s);
		}
		
		reader.close();
	}

	private static Map<String, List<long[]>> parse(File index) throws IOException {
		Map<String, List<long[]>> rtrn=new TreeMap<String, List<long[]>>();
		List<String> lines=BEDFileIO.loadLines(index.getAbsolutePath());
	
		for(String line: lines){
			String[] tokens=line.split("\t");
			String rna=tokens[0];
			List<long[]> list=new ArrayList<long[]>();
			for(int i=1; i<tokens.length; i++){
				long[] array=new long[2];
				array[0]=new Long(tokens[i].split(",")[0]);
				array[1]=new Long(tokens[i].split(",")[1]);
				list.add(array);
			}
			rtrn.put(rna, list);
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		File file=new File("/Users/mguttman/Desktop/head10K");
		String save="/Users/mguttman/Desktop/head10K.index";
		String gene="2010111I01Rik";
		new BuildClusterIndex(file, save);
		readFromIndex(file, new File(save), gene);
		
		
		
		/*String test="Test";
		byte[] b=test.getBytes();
		System.err.println(b.length);
		String test2=new String(b);
		System.err.println(test2);*/
	}
	
}
