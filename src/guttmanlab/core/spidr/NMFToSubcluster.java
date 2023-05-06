package guttmanlab.core.spidr;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class NMFToSubcluster {
	
	private static void writeSubsets(MatrixWithHeaders b, MatrixWithHeaders c, MatrixWithHeaders full, double threshold, String save) throws IOException {
		
		//for each cluster -> get protein list and site list, subset full and write
		for(String cluster: c.getRowNames()) {
			List<String> proteins=getProteins(c, cluster, threshold);
			List<String> sites=getSites(b, cluster, threshold);
			MatrixWithHeaders subset=subset(full, proteins, sites);
			subset.write(save+cluster+".matrix");
		}
		
	}

	

	private static List<String> getProteins(MatrixWithHeaders c, String cluster, double threshold) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String protein: c.getColumnNames()) {
			double score=c.get(cluster, protein);
			if(score>threshold) {rtrn.add(protein);}
		}
		
		return rtrn;
	}



	private static List<String> getSites(MatrixWithHeaders b, String cluster, double threshold) {
		List<String> rtrn=new ArrayList<String>();
		
		for(String site: b.getRowNames()) {
			double score=b.get(site, cluster);
			if(score>threshold) {rtrn.add(site);}
		}
		
		return rtrn;
	}



	private static MatrixWithHeaders subset(MatrixWithHeaders full, List<String> proteins, List<String> sites) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(sites, proteins);
		
		for(String site: sites) {
			for(String protein: proteins) {
				double score=full.get(site, protein);
				rtrn.set(site, protein, score);
			}
		}
		
		return rtrn;
	}



	public static void main(String[] args) throws IOException {
		if(args.length>4) {
		MatrixWithHeaders b=new MatrixWithHeaders(new File(args[0]), ",");
		MatrixWithHeaders c=new MatrixWithHeaders(new File(args[1]), ",");
		MatrixWithHeaders full=new MatrixWithHeaders(new File(args[2]));
		double threshold=Double.parseDouble(args[3]);
		String save=args[4];
		
		writeSubsets(b,c,full, threshold, save);
		}
		else {System.err.println(usage);}
	}

	
	static String usage=" args[0]=basis \n args[1]=coefficient \n args[2]=full matrix \n args[3]=threshold \n args[4]=save";
	
}
