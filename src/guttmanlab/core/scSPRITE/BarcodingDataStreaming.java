package guttmanlab.core.scSPRITE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeSet;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.CombinatoricsUtils;

import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.rnasprite.Cluster;
import guttmanlab.core.simulation.CoordinateSpace;
import htsjdk.samtools.util.CloseableIterator;

public class BarcodingDataStreaming{

		private File barcodeFile;
		private BufferedReader reader;
		private String nextLine;
		Collection<Cluster> clusters;
		
		public BarcodingDataStreaming(File barcodeFile, SingleInterval region) throws IOException{
			this.barcodeFile=barcodeFile;
			this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
			clusters=new ArrayList<Cluster>();
			
			while(hasNext()){
				Cluster c=next();
				if(c.containsOverlappingDNA(region)){clusters.add(c);}
			}
			
		}
		
		public BarcodingDataStreaming(File barcodeFile) throws IOException{
			this.barcodeFile=barcodeFile;
			this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
		}
		
		
		public Collection<Cluster> getClusters(){return this.clusters;}

		
		public boolean hasNext() {
			
			try {
				nextLine = reader.readLine();
					
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				
			}
			return nextLine!=null;
		}


		
		public Cluster next() {
			return Cluster.parseDNACluster(nextLine);
		}

		
		public void close() {
			try {
				reader.close();
				this.reader=new BufferedReader(new InputStreamReader(new FileInputStream(barcodeFile)));
			} catch (IOException e) {
				e.printStackTrace();
			}	
		}
		
		public void reset(){
			close();
		}



		public File getBarcodeFile() {
			return this.barcodeFile;
		}

		public long getNumberOfContacts(int maxClusterSize, int binSize) {
			long sum=0;
			while(hasNext()){
				Cluster c=next();
				Cluster binned=c.bin(binSize);
				int clusterSize=binned.getClusterSize();
				if(clusterSize>1 && clusterSize<maxClusterSize){
					sum+=CombinatoricsUtils.binomialCoefficient(clusterSize, 2);
					//sum+=(clusterSize*clusterSize);
				}
			}
			close();
			return sum;
		}

		public int getNumberOfReads(int maxClusterSize, int binSize) {
			int sum=0;
			while(hasNext()){
				Cluster c=next();
				Cluster binned=c.bin(binSize);
				int clusterSize=binned.getClusterSize();
				if(clusterSize<maxClusterSize){
					sum+=(clusterSize);
				}
			}
			close();
			return sum;
		}

		public Collection<SingleInterval> getFractionsOfBins(int maxClusterSize, int binSize) {
			Collection<SingleInterval> bins=new TreeSet<SingleInterval>();
			//double total=CoordinateSpace.MM10.getBins(binSize).size();
			while(hasNext()){
				Cluster c=next();
				Cluster binned=c.bin(binSize);
				int clusterSize=binned.getClusterSize();
				if(clusterSize<maxClusterSize){
					bins.addAll(binned.getAllDNAIntervals());
				}
			}
			close();
			//double numerator=bins.size();
			//double rtrn=numerator/total;
			return bins;
		}
	
	
}
