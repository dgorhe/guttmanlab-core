package guttmanlab.core.sars.batch;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class BatchGetIntrons {

	
	public static void main(String[] args) throws IOException{
		int size=1000;
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(args[0]);
		String genome=args[1];
		
		
		
		System.err.println("Num introns "+genes.size());
		
		for(int i=0; i<genes.size(); i+=size){
			Collection<Gene> subset=get(genes, i, i+size);
			write(args[2]+".genes"+i+".bed", subset);
			String name="genes"+i;
			getSequences(args[2]+".genes"+i+".bed", args[2]+".genes"+i+".fa", name, genome);
		}
		
	}

	private static Collection<Gene> get(Collection<Gene> genes, int i, int j) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		int count=0;
		System.err.println(i+" "+j+" "+genes.size());
		for(Gene region: genes){
			if(count>=i && count<j){rtrn.add(region);}
			count++;
		}
		System.err.println(rtrn.size());
		return rtrn;
	}

	private static void write(String string, Collection<Gene> subset) throws IOException {
		FileWriter writer =new FileWriter(string);
		
		for(Gene reg: subset){writer.write(reg.toBED()+"\n");}
		
		writer.close();
		
	}

	private static void getSequences(String input, String output, String name, String genome) throws IOException {
		String cmd="java -jar -Xmx8000m /groups/guttman/mguttman/scripts/sars/GetGeneSequences.jar "+input+" "+genome+" "+output;
					
		cmd="srun --nodes=1 --time=24:00:00 --mem=16GB --job-name="+name+ " --output stdout"+name+" --error sterr"+name+ " "+cmd;
			
		Runtime.getRuntime().exec(cmd);
			
		System.err.println(cmd);
		
	}
	
}
