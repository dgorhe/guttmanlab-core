package guttmanlab.core.chipdip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class ProcessSeqFISH {

	public ProcessSeqFISH(File input, Map<String, String> positions, String saveDir) throws IOException {
		Map<String, List<String>> linesByCell=splitByCell(input);
		
		for(String cell: linesByCell.keySet()) {
			String save=saveDir+"/"+cell+".coordinates";
			write(save, linesByCell.get(cell), positions);
			
		}
		
	}

	private Map<String, List<String>> splitByCell(File input) throws IOException {
		Map<String, List<String>> rtrn=new TreeMap<String, List<String>>();
		
		List<String> lines=BEDFileIO.loadLines(input.getAbsolutePath(),1);
		
		for(String line: lines) {
			String[] tokens=line.split(",");
			String cell=tokens[0]+"_"+tokens[1];
			if(!rtrn.containsKey(cell)) {rtrn.put(cell, new ArrayList<String>());}
			List<String> list=rtrn.get(cell);
			list.add(line);
		}
		
		
		return rtrn;
	}

	private void write(String save, List<String> list, Map<String, String> positions) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		List<String> Xs=new ArrayList<String>();
		List<String> Ys=new ArrayList<String>();
		
		for(String line: list) {
			String[] tokens=line.split(",");
			String X=tokens[4];
			String Y=tokens[5];
			Xs.add(X);
			Ys.add(Y);
		}
		
		MatrixWithHeaders matrix=new MatrixWithHeaders(Xs, Ys);
		
		for(String line: list) {
			String[] tokens=line.split(",");
			String name=tokens[3];
			String X=tokens[4];
			String Y=tokens[5];
			String Z=tokens[6];
			int seed=Integer.parseInt(tokens[7]);
			String genomePos=positions.get(name);
			int chrNum=getChrNum(genomePos);
			
			if(seed>=4) {
				writer.write(genomePos+"\t"+name+"\t"+X+"\t"+Y+"\t"+Z+"\n");
				matrix.set(X, Y, chrNum);
			}
			
		}
		
		matrix.write(save+".matrix");
		writer.close();
	}

	
	private int getChrNum(String genomePos) {
		String chrNum=genomePos.split(":")[0].replace("chr", "");
		if(chrNum.equals("X")) {chrNum="20";}
		return Integer.parseInt(chrNum);
	}

	private static Map<String, String> parse(String string) throws IOException {
		List<String> lines=BEDFileIO.loadLines(string,1);
		
		Map<String, String> rtrn=new TreeMap<String, String>();
		for(String line: lines) {
			String[] tokens=line.split("\t");
			String probeName=tokens[1];
			String pos=tokens[3]+":"+tokens[4]+"_"+tokens[5];
			rtrn.put(probeName, pos);
		}
		return rtrn;
		
	}

	
	public static void main(String[] args) throws IOException {
		File input=new File(args[0]);
		Map<String, String> positions=parse(args[1]);
		String saveDir=args[2];
		new ProcessSeqFISH(input, positions, saveDir);
	}

	
	
}
