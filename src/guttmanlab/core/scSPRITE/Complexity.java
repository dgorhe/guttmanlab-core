package guttmanlab.core.scSPRITE;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.simulation.CoordinateSpace;

public class Complexity {
	
	int maxClusterSize=10000;
	public Complexity(File[] files, String save, int binSize) throws IOException{
		int total=CoordinateSpace.MM10.getBins(binSize).size();
		long totalReads=0;
		long totalContacts=0;
		Collection<SingleInterval> allBins=new TreeSet<SingleInterval>();
		FileWriter writer=new FileWriter(save);
		for(int i=0; i<files.length; i++){
			BarcodingDataStreaming data=new BarcodingDataStreaming(files[i]);
			long numContacts=data.getNumberOfContacts(maxClusterSize, binSize);
			int numReads=data.getNumberOfReads(maxClusterSize, binSize);
			totalContacts+=numContacts;
			totalReads+=numReads;
			Collection<SingleInterval> fractionOfBins=data.getFractionsOfBins(maxClusterSize, binSize);
			writer.write(files[i].getName()+"\t"+numReads+"\t"+numContacts+"\t"+fractionOfBins.size()+"\t"+total+"\n");
			if(i%10==0){System.err.println(i);}
			allBins.addAll(fractionOfBins);
		}
		writer.write("all\t"+totalReads+"\t"+totalContacts+"\t"+allBins.size()+"\t"+total);
		writer.close();
	}

	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		int binResolution=new Integer(args[2]);
		new Complexity(files, save, binResolution);
	}
	
}
