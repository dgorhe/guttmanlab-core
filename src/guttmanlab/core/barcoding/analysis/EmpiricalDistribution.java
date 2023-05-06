package guttmanlab.core.barcoding.analysis;

import java.util.ArrayList;
import java.util.Collection;

import guttmanlab.core.annotation.Score;

public class EmpiricalDistribution<T extends Score> {
	
	Collection<T> scores;
	
	public EmpiricalDistribution(){
		this.scores=new ArrayList<T>();
	}

	public void addScore(T score){
		this.scores.add(score);
	}
	
	public double getCDF(T score){
		int count=0;
		int total=0;
		//go through all scores and test those that are less than or equal to score
		for(T val: scores){
			if(val.getScore()<=score.getScore()){count++;}
			total++;
		}
		return ((double)count/(double)total);
	}
	
}
