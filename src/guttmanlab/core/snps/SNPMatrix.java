package guttmanlab.core.snps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.datastructures.Pair;
import guttmanlab.core.math.Statistics;
import jsc.distributions.Hypergeometric;

public class SNPMatrix {

	
	private static File[] getBamFiles(String path) {
		File[] allFiles=new File(path).listFiles();
		
		List<File> files=new ArrayList<File>();
		
		for(int i=0; i< allFiles.length; i++){
			File file=allFiles[i];
			if(file.getName().endsWith(".bam")){files.add(file);}
		}
		
		File[] rtrn=new File[files.size()];
		for(int i=0; i<files.size(); i++){
			rtrn[i]=files.get(i);
		}
		
		return rtrn;
	}
	
	private static MatrixWithHeaders makeMatrix(Map<SingleInterval, Pair<Integer>>[] maps, File[] files) {
		List<String> rows=new ArrayList<String>();
		List<String> columns=new ArrayList<String>();
		
		Collection<SingleInterval> set=new TreeSet<SingleInterval>();
		
		for(int i=0; i< files.length; i++){
			columns.add(files[i].getName()+"_ratio");
			columns.add(files[i].getName()+"_sum");
			for(SingleInterval region: maps[i].keySet()){
				set.add(region);
			}
		}
		
		
		for(SingleInterval region: set){rows.add(region.toUCSC());}
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rows, columns);
		
		for(int i=0; i< maps.length; i++){
			String column1=files[i].getName()+"_ratio";
			String column2=files[i].getName()+"_sum";
			for(SingleInterval region: maps[i].keySet()){
				Pair<Integer> vals=maps[i].get(region);
				String row=region.toUCSC();
				double ratio=getRatio(vals);
				double total=vals.getValue1()+vals.getValue2();
				rtrn.set(row, column1, ratio);
				rtrn.set(row, column2, total);
			}
		}
		
		return rtrn;
	}
	
	private static double getRatio(Pair<Integer> vals) {
		double sum=vals.getValue1()+vals.getValue2();
		return (double)vals.getValue1()/sum;
	}
	
	private static MatrixWithHeaders filter(MatrixWithHeaders mwh, int minCount) {
		Collection<String> columns=new TreeSet<String>();
		for(String column: mwh.getColumnNames()){
			if(column.endsWith("sum")){
				columns.add(column);
			}
		}
		
		MatrixWithHeaders temp=mwh.submatrixByColumnNames(columns);
		
		Collection<String> rows=new TreeSet<String>();
		
		for(String row: temp.getRowNames()){
			double[] vals=temp.getRow(row);
			double max=Statistics.max(vals);
			if(max>minCount){rows.add(row);}
		}
		
		return mwh.submatrixByRowNames(rows);
		
	}
	
	
	
	private static MatrixWithHeaders computeHypergeometric(Map<SingleInterval, Pair<Integer>>[] maps, Map<SingleInterval, Pair<Integer>> controlScores, MatrixWithHeaders mwh, int minCount, String controlName, File[] files) {
		MatrixWithHeaders rtrn=initialize(mwh.getRowNames(), getNames(files), 1.0);
		
		Collection<String> rows=new TreeSet<String>();
		for(int i=0; i<files.length; i++){
			String sample=files[i].getName()+"_p";
			String sampleFold=files[i].getName()+"_fold";
			
			Map<SingleInterval, Pair<Integer>> sampleScores=maps[i];
			for(SingleInterval region: sampleScores.keySet()){
				if(rtrn.containsRow(region.toUCSC())){
					Pair<Integer> sampleVal=sampleScores.get(region);
					Pair<Integer> nonSampleVals=get(maps, i, region);
					if(controlScores.containsKey(region)){
						Pair<Integer> controlVal=controlScores.get(region);
						double p=hypergeometric(sampleVal, controlVal, minCount);
						double p2=hypergeometric(sampleVal, nonSampleVals, minCount);
						p=Math.max(p, p2);
						
						double expected=((double)nonSampleVals.getValue1()/(sum(nonSampleVals)))*sum(sampleVal);
						double enrich=(double)sampleVal.getValue1()/expected;
						
						rtrn.set(region.toUCSC(), sample, p);
						rtrn.set(region.toUCSC(), sampleFold, enrich);
						rows.add(region.toUCSC());
					}
				}
			}
			
		}
		
		rtrn=rtrn.submatrixByRowNames(rows);
		//rtrn=filterMinPValue(rtrn);
		return rtrn;
		
	}

	private static Pair<Integer> get(Map<SingleInterval, Pair<Integer>>[] maps, int pos, SingleInterval region) {
		Pair<Integer> rtrn=new Pair<Integer>(0,0);
		
		for(int i=0; i<maps.length; i++){
			if(i!=pos){
				if(maps[i].containsKey(region)){
					Pair<Integer> val=maps[i].get(region);
					int val1=rtrn.getValue1()+val.getValue1();
					int val2=rtrn.getValue2()+val.getValue2();
					rtrn.setValue1(val1);
					rtrn.setValue2(val2);
				}
			}
		}
		
		return rtrn;
	}

	private static MatrixWithHeaders filterMinPValue(MatrixWithHeaders rtrn) {
		Collection<String> names=new TreeSet<String>();
		for(String row: rtrn.getRowNames()){
			double[] vals=rtrn.getRow(row);
			double min=Statistics.min(vals);
			if(min<=0.95){names.add(row);}
		}
		return rtrn.submatrixByRowNames(names);
	}

	private static MatrixWithHeaders initialize(List<String> rowNames, List<String> names, double d) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rowNames, names);
		
		
		for(String row: rtrn.getRowNames()){
			for(String column: rtrn.getColumnNames()){
				rtrn.set(row, column, d);
			}
		}
		
		return rtrn;
	}

	private static double hypergeometric(Pair<Integer> sampleVal, Pair<Integer> controlVal, int minCount) {
		ChiSquareTest test=new ChiSquareTest();
		
		long[] observed1=new long[2];
		long[] observed2=new long[2];
		
		observed1[0]=sampleVal.getValue1();
		observed1[1]=sampleVal.getValue2();
		
		observed2[0]=controlVal.getValue1();
		observed2[1]=controlVal.getValue2();
		
		double ratio1=getRatio(sampleVal);
		double ratio2=getRatio(controlVal);
		
		double p=1.0;
		if(ratio1==0 && ratio2==0){
			p=1.0;
		}
		else if(ratio1==1 && ratio2==1){
			p=1.0;
		}
		else if(sum(sampleVal)>minCount && sum(controlVal)>minCount){
			//System.err.println(observed1[0]+" "+observed1[1]+" "+observed2[0]+" "+observed2[1]);
			p=test.chiSquareTestDataSetsComparison(observed1, observed2);
		}
		
		
		return p;
	}

	private static int sum(Pair<Integer> observed1) {
		return observed1.getValue1()+observed1.getValue2();
	}

	private static List<String> getNames(File[] files) {
		List<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<files.length ; i++){
			rtrn.add(files[i].getName()+"_p");
			rtrn.add(files[i].getName()+"_fold");
		}
		
		return rtrn;
	}

	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File[] files=getBamFiles(args[0]);
			String saveDir=args[1];
			File snpFile=new File(args[2]);
			int minCount=new Integer(args[3]);
			File control=new File(args[4]);
			
			Map<String, IntervalTree<String>> snpTree=QuantifyBySNP.parseTree(snpFile);
			
			Map<SingleInterval, Pair<Integer>>[] maps=new Map[files.length];
			
			for(int i=0; i<files.length; i++){
				File bamFile=files[i];
				System.err.println(bamFile.getAbsolutePath());
				String save=saveDir+"/"+files[i].getName();
				QuantifyBySNP snp=new QuantifyBySNP(bamFile, snpTree, save, true);
				Map<SingleInterval, Pair<Integer>> snpCounts=snp.getSNPCounts();
				maps[i]=snpCounts;
			}
			
			QuantifyBySNP snp=new QuantifyBySNP(control, snpTree, saveDir+"/"+control.getName(), true);
			Map<SingleInterval, Pair<Integer>> sampleMap=snp.getSNPCounts();
			
			MatrixWithHeaders mwh=makeMatrix(maps, files);
			mwh=filter(mwh, minCount);
			
			MatrixWithHeaders hyper=computeHypergeometric(maps, sampleMap, mwh, minCount, control.getName(), files);
			
			
			mwh.write(saveDir+"/snps.counts");
			hyper.write(saveDir+"/snps.hypergeometric");
		}
		else{System.err.println(usage);}
	}

	

	



	static String usage=" args[0]=folder of bams \n args[1]=save directory \n args[2]=snp file \n args[3]=min count \n args[4]=control file";
	
	
	
}
