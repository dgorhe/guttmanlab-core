package guttmanlab.core.sars.leader;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.math3.stat.inference.ChiSquareTest;

import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;

public class VariantsPerVirus {
	
	Collection<String> proteins;
	
	public VariantsPerVirus (String file, String save, int cutoff) throws IOException {
		this.proteins=new TreeSet<String>();
		this.proteins.add("NSP1");
		this.proteins.add("Spike");
		
		//Go through each virus (line) and get all variants
		//Make matrix where rows are viruses and columns positions in protein
		List<String> lines=BEDFileIO.loadLines(file);
		
		Map<String, Double> freq=getFrequencies(lines);
		
		System.err.println("loaded");
		List<String> positions=getPositions(lines);
		MatrixWithHeaders observed=scoreObserved(lines, positions, freq);
		observed.write(save);
		
		observed=filter(observed, cutoff, (0.1*cutoff));
		observed.write(save+".filtered");
		
		MatrixWithHeaders pvals=computePVals(observed, freq, lines.size()-1);
		pvals.write(save+".pvals");
		
		/*MatrixWithHeaders norm=normalize(observed, freq, lines.size()-1);
		norm.write(save+".normalized");
		norm=filter(norm);*/
	}
	
	
	private MatrixWithHeaders computePVals(MatrixWithHeaders observed, Map<String, Double> freq, int total) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		for(String row: observed.getRowNames()) {
			for(String column: observed.getColumnNames()) {
				if(!row.equals(column)) {
					double ab=observed.get(row, column);
					double a=freq.get(row);
					double b=freq.get(column);
					double aNb=a-ab;
					double bNa=b-ab;
					double nbNa=total-a-b+ab;
					
					//double p=Statistics.fisherExact(ab, aNb, bNa, nbNa);
					
					
					
					long[][] matrix=new long[2][2];
					matrix[0][0]=(long)ab;
					matrix[0][1]=(long)aNb;
					matrix[1][0]=(long)bNa;
					matrix[1][1]=(long)nbNa;
					
					ChiSquareTest chi=new ChiSquareTest();
					double p=chi.chiSquareTest(matrix);
					
					
					double freqAB=ab/(double)total;
					double freqA=a/(double)total;
					double freqB=b/(double)total;
					double expected=freqA*freqB*(double)total;
					double enrichment=Math.log(freqAB/(freqA*freqB))/Math.log(2);
					
					System.err.println(row+" "+column+" "+ab+" "+ expected+" "+a+" "+b+" "+enrichment+" "+p);
					
					//double score=-Math.log(p)/Math.log(10);
					//if(ab<expected) {score=-score;}
					if(p<0.001 && Math.max(ab, expected)>10) {
						rtrn.set(row,column, enrichment);
					}
				}
				
			}
		}
		System.err.println("done");
		
		return rtrn;
	}


	private MatrixWithHeaders filter(MatrixWithHeaders observed, int cutoff, double min) {
		Collection<String> list=new TreeSet<String>();
		for(String m: observed.getRowNames()) {
			double val=observed.get(m, m);
			double[] vals=observed.getColumn(m);
			double max=Statistics.max(vals, val);
			if(val>cutoff && max>min) {list.add(m);}
		}
		
		observed=observed.submatrixByColumnNames(list);
		observed=observed.submatrixByRowNames(list);
		
		return observed;
	}


	private MatrixWithHeaders normalize(MatrixWithHeaders observed, Map<String, Double> freqs, int total) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(observed.getRowNames(), observed.getColumnNames());
		
		
		for(String row: observed.getRowNames()) {
			for(String column: observed.getColumnNames()) {
				double o=observed.get(row, column);
				if(o>10) {
					double freqAB=o/(double)total;
					double freqA=freqs.get(row)/(double)total;
					double freqB=freqs.get(column)/(double)total;
					double enrichment=Math.log(freqAB/(freqA*freqB))/Math.log(2);
					System.err.println(row+" "+column+" "+o+" "+freqAB+" "+enrichment);
					rtrn.set(row, column, enrichment);
				}
			}
		}
		
		
		return rtrn;
	}

	
	private Map<String, Double> getFrequencies(List<String> lines) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		
		int counter=0;
		for(String line: lines) {
			Collection<String> mutants=getMutants(line);
			for(String m1: mutants) {
				double count=0;
				if(rtrn.containsKey(m1)) {count=rtrn.get(m1);}
				count++;
				rtrn.put(m1, count);
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		
		return rtrn;
	}
	

	private MatrixWithHeaders scoreObserved(List<String> lines, List<String> positions, Map<String, Double> freq) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(positions, positions);
		
		
		int counter=0;
		for(String line: lines) {
			Collection<String> mutants=getMutants(line);
			for(String m1: mutants) {
				for(String m2: mutants) {
					if(!m1.equals(m2) && positions.contains(m1) && positions.contains(m2)) {
						rtrn.incrementCount(m1, m2);
					}
				}
			}
			counter++;
			if(counter%100000==0) {System.err.println(counter);}
		}
		
		
		for(String m1: positions) {
			rtrn.set(m1, m1, freq.get(m1));
		}
		
		
		return rtrn;
	}


	private MatrixWithHeaders filter(MatrixWithHeaders mwh) {
		List<String> rows=new ArrayList<String>();
		for(String row: mwh.getRowNames()) {
			double[] vals=mwh.getRow(row);
			double max=Statistics.maxAbs(vals);
			if(max>3) {rows.add(row);}
		}
		
		mwh=mwh.submatrixByRowNames(rows);
		mwh=mwh.submatrixByColumnNames(rows);
		
		return mwh;
	}


	private Collection<String> getMutants(String line) {
		String[] tokens=line.split("\t");
		Collection<String> mutants=parseMutants(tokens[5]);
		return mutants;
	}




	
	
	private List<String> getPositions(List<String> lines) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String line: lines) {
			String[] tokens=line.split("\t");
			
			
			Collection<String> mutants=parseMutants(tokens[5]);
			rtrn.addAll(mutants);
			
		}
		
		List<String> list=new ArrayList<String>();
		
		for(String m: rtrn) {
			String protein=m.split("_")[0];
			if(this.proteins.contains(protein)) {
				list.add(m);
			}
		}
		
		
		return list;		
	}
	


	private Collection<String> parseMutants(String string) {
		Collection<String> rtrn=new TreeSet<String>();
		
		if(string.equals("()") || string.isEmpty()){return rtrn;}
		else if(!string.contains("(")) {return rtrn;}
		
		String trimmed=string.split("\\(")[1].split("\\)")[0];
			
		
		
		
		String[] tokens=trimmed.split(",");
		
		for(int i=0; i<tokens.length; i++) {
			
			if(tokens[i].split("_").length<2) {System.err.println(tokens[i]);} //TODO what's going on?
			else {
				String protein=tokens[i].split("_")[0];
				String pos=getNum(tokens[i].split("_")[1]);
				rtrn.add(protein+"_"+pos);
			}
		}
		
		return rtrn;
	}
	



	private String getNum(String string) {
		String rtrn=string;
		if(string.contains("del")) {
			rtrn=string.substring(1, string.length()-3);
		}
		else if(string.contains("stop")) {
			rtrn=string.substring(1, string.length()-4);
		}
		else {
			rtrn=string.substring(1, string.length()-1);
			//System.err.println(string+" "+rtrn);
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws IOException {
		new VariantsPerVirus(args[0], args[1], Integer.parseInt(args[2]));
	}
	
}
