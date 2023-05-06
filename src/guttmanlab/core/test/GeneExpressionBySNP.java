package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.util.GeneExpressionUtils;

public class GeneExpressionBySNP {

	
	
	public static void main(String[] args) throws IOException{
		File[] bams=getBAMs(new File(args[0]).listFiles());
		String genes=(args[1]);
		String save=args[2];
		
		
		MatrixWithHeaders mwh=GeneExpressionUtils.quantify(bams, genes);
		
		mwh=filter(mwh);
		
		mwh.write(save);
	}

	private static MatrixWithHeaders filter(MatrixWithHeaders mwh) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()){
			double[] vals=mwh.getRow(row);
			if(Statistics.sum(vals)>0){list.add(row);}
		}
		return mwh.submatrixByRowNames(list);
	}

	private static File[] getBAMs(File[] listFiles) {
		int count=0;
		for(int i=0; i<listFiles.length; i++){
			if(listFiles[i].getName().endsWith("bam")){count++;}
		}
		
		File[] rtrn=new File[count];
		
		int counter=0;
		for(int i=0; i<listFiles.length; i++){
			if(listFiles[i].getName().endsWith("bam")){
				rtrn[counter]=listFiles[i];
				counter++;
			}
		}
		
		return rtrn;
	}
	
}
