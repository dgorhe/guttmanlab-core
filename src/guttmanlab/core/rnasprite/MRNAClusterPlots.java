package guttmanlab.core.rnasprite;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class MRNAClusterPlots {

	public MRNAClusterPlots(BarcodingDataStreaming allData, String rna, Collection<String> mRNAs, List<String> geneSet, String save) throws IOException{
		Collection<Cluster> clusters=allData.getRNAClusters(rna);
		
		System.err.println("Got all clusters");
		
		Collection<Cluster> hasMRNA=hasMRNA(clusters, mRNAs);
		MatrixWithHeaders withMRNA=allData.getClusterRNAInteractionMatrix(hasMRNA, geneSet);
		withMRNA.write(save+".withMRNA");
		
		System.err.println("with RNA made");
		
		Collection<Cluster> noMRNA=noMRNA(clusters, mRNAs);
		MatrixWithHeaders withoutMRNA=allData.getClusterRNAInteractionMatrix(noMRNA, geneSet);
		withoutMRNA.write(save+".withoutMRNA");
		
		System.err.println("without RNA made");
		
	}
	
	private Collection<Cluster> noMRNA(Collection<Cluster> clusters, Collection<String> mRNAs) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		for(Cluster c: clusters){
			boolean hasMRNA=hasMRNA(c, mRNAs);
			if(!hasMRNA){rtrn.add(c);}
		}
		
		return rtrn;
	}



	private Collection<Cluster> hasMRNA(Collection<Cluster> clusters, Collection<String> mRNAs) {
		Collection<Cluster> rtrn=new ArrayList<Cluster>();
		
		for(Cluster c: clusters){
			boolean hasMRNA=hasMRNA(c, mRNAs);
			if(hasMRNA){rtrn.add(c);}
		}
		
		return rtrn;
	}


	private static List<String> parse(String string) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		List<String> lines=BEDFileIO.loadLines(string);
		
		for(String line: lines){
			rtrn.add(line);
		}
		
		return rtrn;
	}

	private boolean hasMRNA(Cluster c, Collection<String> mRNAs) {
		for(String gene: c.getRNANames()){
			if(mRNAs.contains(gene)){return true;}
		}
		return false;
	}

	public static void main(String[] args) throws IOException{
		BarcodingDataStreaming data=new BarcodingDataStreaming(new File(args[0]));
		String rna=args[1];
		Collection<String> mRNAs=parse(args[2]);
		List<String> genesToUse=parse(args[3]);
		String save=args[4];
		new MRNAClusterPlots(data, rna, mRNAs, genesToUse, save);
	}
	
	
}
