package guttmanlab.core.primer3;


import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collection;

import org.apache.log4j.Logger;


import jaligner.Alignment;

public class PrimerUtils {
	
	private static Logger logger = Logger.getLogger(PrimerUtils.class.getName());
	private static final float MAX_PRIMER_PENALTY = (float)0.05;

	
	/**
	 * Get one primer pair
	 * @param primerLength Primer length
	 * @param pathPrimer3core primer3core executable
	 * @param existingPrimerReader Reader for a table of existing primer pairs, each line having fields as defined in the constructor PrimerPair(String[]), or null if not using
	 * @param primerId Primer ID or null
	 * @return Primer pair with less than max primer penalty
	 * @throws IOException
	 */
	/*public static PrimerPair getOneSyntheticPrimerPair(int primerLength, String pathPrimer3core, double optimalMeltingTemp, BufferedReader existingPrimerReader, String primerId) throws IOException {
		if(existingPrimerReader == null) {
			return getOneSyntheticPrimerPair(primerLength, pathPrimer3core, optimalMeltingTemp, primerId);
		}
		StringParser s = new StringParser();
		while(existingPrimerReader.ready()) {
			s.parse(existingPrimerReader.readLine());
			PrimerPair p = new PrimerPair(s.getStringArray());
			// Check that the primer pair read from the file satisfies the constraints
			if(p.getLeftPrimer().length() != primerLength || p.getRightPrimer().length() != primerLength) {
				logger.warn("Skipping existing primer pair " + p.getLeftPrimer() + " " + p.getRightPrimer() + " because requested length is " + primerLength);
				continue;
			}
			if(p.getLeftPrimerTM() < optimalMeltingTemp - 1 || p.getRightPrimerTM() > optimalMeltingTemp + 1) {
				logger.warn("Skipping existing primer pair " + p.getLeftPrimer() + " " + p.getRightPrimer() + " because left primer TM is " + p.getLeftPrimerTM());
				continue;
			}
			if(p.getRightPrimerTM() < optimalMeltingTemp - 1 || p.getRightPrimerTM() > optimalMeltingTemp + 1) {
				logger.warn("Skipping existing primer pair " + p.getLeftPrimer() + " " + p.getRightPrimer() + " because right primer TM is " + p.getRightPrimerTM());
				continue;
			}
			if(primerId != null) p.setId(primerId);
			return p;
		}
		//logger.warn("Ran out of existing primer pairs in file. Creating new primer pair.");
		return getOneSyntheticPrimerPair(primerLength, pathPrimer3core, optimalMeltingTemp, primerId);
	}*/
	
	/**
	 * Get one primer pair
	 * @param primerLength Primer length
	 * @param pathPrimer3core primer3core executable
	 * @param primerId Primer ID or null
	 * @return Primer pair with less than max primer penalty
	 * @throws IOException
	 */
	/*public static PrimerPair getOneSyntheticPrimerPair(int primerLength, String pathPrimer3core, double optimalMeltingTemp, String primerId) throws IOException {
		int numTried = 0;
		// Repeat until a suitable primer pair is found
		while(true) {
			numTried++;
			String seq = Sequence.generateRandomSequence(5000);
			// Only ask for one primer pair
			Collection<PrimerPair> primers = PcrPrimerDesigner.designSyntheticPrimers(seq, 1, primerLength, 5000, pathPrimer3core, optimalMeltingTemp);
			// A primer pair was found
			if(primers != null && !primers.isEmpty()) {
				PrimerPair primer = primers.iterator().next();
				if(primer.getPrimerPairPenalty() <= MAX_PRIMER_PENALTY) {
					if(primerId != null) primer.setId(primerId);
					return primer;
				}
			}
			if(numTried == 1000) {
				throw new IllegalStateException("Tried 1000 random primer pairs without finding an acceptable one. Try adjusting primer length and optimal Tm.");
			}
			if(numTried % 100 == 0) {
				logger.info("Tried " + numTried + " random primer pairs");
			}
		}
	}*/

	
	/*public static double computeMaxTM(Alignment align){
		return Math.max(computeTM(align), computeTM2(align));
	}*/

	public static double computeTM(String seq){
		//Tm = 64.9C + 41C x (number of Gs and Cs in the primer  16.4)/N
		char[] bases=seq.toCharArray();
		double TM=64.9;
		int numGC=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='G'|| bases[i]=='C' || bases[i]=='c' || bases[i]=='g'){numGC++;}
		}
		TM=TM+41*((numGC-16.4)/(double)bases.length);
		
		return TM;
	}
	
	
	public static double percentGC(String seq) {
		char[] bases=seq.toCharArray();
		int numGC=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='G'|| bases[i]=='C' || bases[i]=='c' || bases[i]=='g'){numGC++;}
		}
				
		return 100.0*(numGC/(double)bases.length);
	}

	public static double computeTM(Alignment align) {
		//COmpute the TM only across the matches
		String concat="";
		//Tm = 64.9C + 41C x (number of Gs and Cs in the primer  16.4)/N
		double TM=64.9;
		int numGC=0;
		int numMismatch=0;
		double length=0.0;
		for(int i=0; i<align.getSequence1().length; i++){
			char base1=align.getSequence1()[i];
			char base2=align.getSequence2()[i];
			if(base1==base2 && base1!='-'){
				if(base1=='G'|| base1=='C' || base1=='c' || base1=='g'){numGC++;}
				concat=concat+base1;
				length++;
			}
			else if(base1!=base2 && base1!='-' && base2!='-'){
				numMismatch++;
			}
		}
		
		TM=TM+41*((numGC-16.4)/length);
		double percentMismatch=((numMismatch/(length+numMismatch)))*100;
		double TM1=Math.max(0, TM-percentMismatch);
		//System.err.println(TM+" "+TM1+" "+numMismatch+" "+length+" "+percentMismatch);
		return TM1;
	}
	
	public static void main(String[] args){
		String seq="TAATACGACTCACTATAGGG";
		PrimerPair p=new PrimerPair();
		
		System.err.println(seq+" "+computeTM(seq));
	}


}
