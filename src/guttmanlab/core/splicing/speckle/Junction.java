package guttmanlab.core.splicing.speckle;

import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;

public class Junction implements Comparable<Junction>{
	SingleInterval exon1;
	SingleInterval exon2;
	SingleInterval tss;
	SingleInterval intron;

	public Junction(Gene j, int tss2) {
		this.intron=j.getIntrons().iterator().next().getSingleInterval();
		this.exon1=j.getFirstBlock();
		this.exon2=j.getLastBlock();
		this.tss=new SingleInterval(j.getReferenceName(), tss2, tss2+1);
	}

	@Override
	public int compareTo(Junction o) {
		int c1=exon1.compareTo(o.exon1);
		int c2=exon2.compareTo(o.exon2);
		int c3=intron.compareTo(o.intron);
		int c4=tss.compareTo(o.tss);
		
		if(c1!=0) {return c1;}
		if(c2!=0) {return c2;}
		if(c3!=0) {return c3;}
		return c4;
	}
	
	
	
	
	
	
}
