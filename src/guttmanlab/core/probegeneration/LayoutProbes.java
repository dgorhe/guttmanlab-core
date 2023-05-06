package guttmanlab.core.probegeneration;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.Pair;

public class LayoutProbes {

	public LayoutProbes(Collection<Map<String, String>> probes, Collection<String> primers, String save, Collection<String> filterList) throws IOException{
		FileWriter writer=new FileWriter(save+".fa");
		
		Map<String, String> probesWithPrimers=addPrimers(probes, primers, writer, filterList);
		
		writer.close();
		write(save+".probes", probesWithPrimers);
		
		Collection<String> finalProbes=assembleIntoProbes(probesWithPrimers);
		
		write(save+".finalOligos", finalProbes);
	}

	private void write(String save, Map<String, String> probesWithPrimers) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String probe: probesWithPrimers.keySet()){
			writer.write(probe+"\t"+probesWithPrimers.get(probe)+"\n");
		}
		
		writer.close();
		
	}

	private Map<String, String> addPrimers(Collection<Map<String, String>> probes, Collection<String> primers, FileWriter writer, Collection<String> filterList) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		Map<String, Pair<String>> primerSet=getPrimers(primers);
		
		int counter=0;
		for(Map<String, String> set: probes){
			String setName="set"+counter;
			System.err.println(setName);
			Pair<String> pp=primerSet.get(setName);
			for(String probe: set.keySet()){
				String pos=set.get(probe);
				String name=pos+"_"+setName;
				if(!filterList.contains(name)){
					String fullProbe=pp.getValue1()+probe+pp.getValue2();
					rtrn.put(fullProbe, setName);
					writer.write(">"+name+"\n");
					writer.write(fullProbe+"\n");
				}
			}
			counter++;
		}
		return rtrn;
	}

	private Map<String, Pair<String>> getPrimers(Collection<String> primers) {
		Map<String, Pair<String>> rtrn=new TreeMap<String, Pair<String>>();
		
		int counter=0;
		for(String primer: primers){
			String setName="set"+counter;
			Pair<String> pp=new Pair<String>(primer.split("\t")[2], primer.split("\t")[4]);
			rtrn.put(setName, pp);
			counter++;
		}
		
		return rtrn;
	}

	private Collection<String> assembleIntoProbes(Map<String, String> probesWithPrimers) {
		Collection<String> rtrn=new TreeSet<String>();
		
		Collection<String> residualList=new TreeSet<String>();
		residualList.addAll(probesWithPrimers.keySet());
		
		for(String probe: probesWithPrimers.keySet()){
			if(residualList.contains(probe)){
				String set=probesWithPrimers.get(probe);
				residualList.remove(probe);
				String probe2=getNext(residualList, set, probesWithPrimers);
				String mergedProbe=merge(probe, probe2);
				rtrn.add(mergedProbe);
				residualList.remove(probe2);
			}
		}
		
		return rtrn;
	}

	private String getNext(Collection<String> residualList, String set, Map<String, String> probesWithPrimers) {
		for(String probe2: residualList){
			String set2=probesWithPrimers.get(probe2);
			if(!set.equalsIgnoreCase(set2)){return probe2;}
		}
		return "";
	}

	private String merge(String probe, String probe2) {
		return probe+probe2;
	}

	private void write(String save, Collection<String> finalProbes) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String probe: finalProbes){writer.write(probe+"\n");}
		
		writer.close();
	}
	
	private static Map<String, String> parseProbes(File string) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		List<String> list=BEDFileIO.loadLines(string.getAbsolutePath());
		for(String line: list){
			String probe=line.split("\t")[3];
			String pos=line.split("\t")[0]+":"+line.split("\t")[1]+"-"+line.split("\t")[2];
			rtrn.put(probe, pos);
		}
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		Collection<String> primers=BEDFileIO.loadLines(args[1]);
		String save=args[2];
		Collection<String> filterList=parseFilter(args[3]);
		
		Collection<Map<String, String>> allProbes=new ArrayList<Map<String, String>>();
		for(int i=0; i<files.length; i++){
			System.err.println(files[i].getAbsolutePath()+" set "+i);
			Map<String, String> probes=parseProbes(files[i]);
			allProbes.add(probes);
		}
		new LayoutProbes(allProbes, primers, save, filterList);
	}

	private static Collection<String> parseFilter(String string) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		
		List<String> lines=BEDFileIO.loadLines(string);
		for(String line: lines){
			rtrn.add(line.split("\t")[0]);
		}
		
		return rtrn;
	}

	
	
}
