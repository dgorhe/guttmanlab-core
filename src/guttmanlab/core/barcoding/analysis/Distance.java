package guttmanlab.core.barcoding.analysis;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;


public class Distance implements Comparable<Distance>{
		int[] distanceList;
		int resolution;
		int genomeLength;
		
		
		public Distance(int[] distanceList, int resolution){
			Arrays.sort(distanceList);
			this.distanceList=distanceList;
			this.resolution=resolution;
			this.genomeLength=computeGenomeLength();
		}

		
		
		private int computeGenomeLength() {
			int sum=0;
			for(int i=0; i<distanceList.length; i++){
				sum+=distanceList[i];
			}
			int blockCount=(distanceList.length+1)*resolution;
			return sum+blockCount;
		}
		
		public int getTotalGenomeLength() {
			return this.genomeLength;
		}

		public int getResolution() {
			return this.resolution;
		}

		public void setResolution(Integer resolution) {
			this.resolution=resolution;
		}

		@Override
		public int compareTo(Distance o) {
			int sizes=this.resolution-o.resolution;
			if(sizes!=0){return sizes;}
			
			sizes=this.distanceList.length-o.distanceList.length;
			if(sizes != 0){return sizes;}
			
			//go element by element and test if equal
			for(int i=0; i<this.distanceList.length; i++){
				sizes=this.distanceList[i]-o.distanceList[i];
				if(sizes!=0){return sizes;}
			}
			return 0;
		}
		
		@Override
		public String toString(){
			String rtrn="[";
			
			for(int i=0; i<this.distanceList.length; i++){
				rtrn+=this.distanceList[i];
				if(i!=this.distanceList.length-1){rtrn+=",";}
			}
			rtrn+="] ";
			
			return rtrn;	
		}

		public String toFileName() {
			String rtrn="r"+resolution;
			for(int i=0; i<distanceList.length; i++){
				rtrn+="_"+distanceList[i];
			}
			return rtrn;
		}

		public int[] getPermutedDistanceList() {
			ArrayList<Integer> list=new ArrayList<Integer>();
			for(Integer distance: distanceList){list.add(distance);}
			int[] rtrn=new int[distanceList.length];
			for(int i=0; i<rtrn.length; i++){
				int randomIndex=new Double(Math.random()*list.size()).intValue();
				rtrn[i]=list.remove(randomIndex);
			}
			return rtrn;
		}

		
	}
	
