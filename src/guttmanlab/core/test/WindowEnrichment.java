package guttmanlab.core.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.datastructures.MatrixWithHeaders;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.simulation.CoordinateSpace;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class WindowEnrichment {

	MatrixWithHeaders mwh;
	
	public WindowEnrichment(File[] bamFiles, int binResolution) throws IOException{
		List<String> samples=new ArrayList<String>();
		for(int i=0; i<bamFiles.length; i++){
			if(bamFiles[i].getName().endsWith(".bam")){
				samples.add(bamFiles[i].getName());
			}
		}
		
		MatrixWithHeaders mwh=new MatrixWithHeaders(CoordinateSpace.MM10.getBins(binResolution), samples);
		
		for(int i=0; i<bamFiles.length; i++){
			if(bamFiles[i].getName().endsWith(".bam")){
				score(bamFiles[i], mwh, binResolution);
			}
		}
		
		this.mwh=filter(mwh);
		//mwh.write(save);
		
	}

	private MatrixWithHeaders filter(MatrixWithHeaders mwh) {
		List<String> list=new ArrayList<String>();
		for(String row: mwh.getRowNames()){
			double[] vals=mwh.getRow(row);
			if(hasVal(vals)){list.add(row);}
		}
		
		return mwh.submatrixByRowNames(list);
	}

	private boolean hasVal(double[] vals) {
		for(int i=0; i<vals.length; i++){
			if(vals[i]<=0){return false;}
		}
		return true;
	}
	
	
	private static MatrixWithHeaders add(File bamFile, MatrixWithHeaders mwh1, int binResolution) {
		List<String> list=new ArrayList<String>();
		list.add(bamFile.getName());
		
		MatrixWithHeaders temp=new MatrixWithHeaders(mwh1.getRowNames(), list);
		
		
	
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			SingleInterval bin=f.bin(binResolution);
			if(temp.containsRow(bin.toUCSC())){
				temp.incrementCount(bin.toUCSC(), bamFile.getName());
			}
			
			counter++;
			
		}
		
		reads.close();
		reader.close();
		mwh1.appendColumns(temp);
		
		List<String> toUse=new ArrayList<String>();
		for(String row: mwh1.getRowNames()){
			if(Statistics.max(mwh1.getRow(row))>10){toUse.add(row);}
		}
		
		mwh1=mwh1.submatrixByRowNames(toUse);
		
		return mwh1;
	}

	private void score(File bamFile, MatrixWithHeaders mwh, int binResolution) {
		SAMFileReader reader=new SAMFileReader(bamFile);
		SAMRecordIterator reads=reader.iterator();
		
		int counter=0;
		while(reads.hasNext()){
			SAMRecord record=reads.next();
			
			SAMFragment f=new SAMFragment(record);
			SingleInterval bin=f.bin(binResolution);
			if(mwh.containsRow(bin.toUCSC())){
				mwh.incrementCount(bin.toUCSC(), bamFile.getName());
				//System.err.println(counter+" "+bin.toUCSC() +" "+ bamFile.getName()+" "+ mwh.get(bin.toUCSC(), bamFile.getName()));
			}
			
			counter++;
			
		}
		
		reads.close();
		reader.close();
	}
	
	public static void main(String[] args) throws IOException{
		File[] files=new File(args[0]).listFiles();
		int binResolution=new Integer(args[1]);
		String save=args[2];
		WindowEnrichment w=new WindowEnrichment(files, binResolution);
		MatrixWithHeaders mwh=w.mwh;
		mwh=add(new File(args[3]), mwh, binResolution);
		mwh.write(save);
	}

	
	
}
