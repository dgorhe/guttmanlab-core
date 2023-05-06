package guttmanlab.core.rnasprite.hubs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;

public class SubsetHeatmap {

	
	private static void subset(MatrixWithHeaders mwh, Collection<SingleInterval> regions, int binResolution, String save) throws IOException {
		
		List<String> bins=bin(regions, binResolution, mwh);
		
		mwh=mwh.submatrixByRowNames(bins);
		mwh=mwh.submatrixByColumnNames(bins);
		
		mwh.write(save);
		
	}

	private static List<String> bin(Collection<SingleInterval> regions, int binResolution, MatrixWithHeaders mwh) {
		List<String> rtrn=new ArrayList<String>();
		Collection<SingleInterval> temp=new TreeSet<SingleInterval>();
		
		for(SingleInterval r: regions) {
			Collection<SingleInterval> bins=r.allBins(binResolution);
			temp.addAll(bins);
		}
		
		
		for(SingleInterval r: temp) {
			String name=r.toUCSC();
			if(mwh.containsColumn(name)) {
				rtrn.add(name);
			}
			else {System.err.println(name);}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws IOException {
		MatrixWithHeaders mwh=new MatrixWithHeaders(new File(args[0]));
		Collection<SingleInterval> genes=BEDFileIO.loadSingleIntervalFromFile(args[1]);
		int binResolution=1000000;
		String save=args[2];
		subset(mwh, genes, binResolution, save);
	}
	
}
