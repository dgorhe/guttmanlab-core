package guttmanlab.core.barcoding.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;


import org.paukov.combinatorics3.Generator;

import guttmanlab.core.simulation.CoordinateSpace;

public class DistanceScoring {
	Integer[] chromosomalDistancesArray;
	ArrayList<Integer> chromosomalDistancesList;

	public DistanceScoring(int nmerSize, int resolution, String type) {
		//this should probably changed but we're storing chromosome sizes in memory
		
		if (type.equals("H")) {
			// hg19
			Integer[] temp = { 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
					141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210,
					78077248, 63025520, 59128983, 51304566, 48129895, 155270560, 59373566 };
			chromosomalDistancesArray = temp;
		} else if (type.equals("M")) {
			// mm10
			/*Integer[] temp = { 195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213,
					124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271,
					90702639, 61431566, 171031299, 91744698,16299 };*/
			
			
			// mm9
            Integer[] temp = { 197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871,
                    124076172, 129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651,
                    90772031, 61342430, 166650296, 15902555, 16299 };
            
			chromosomalDistancesArray = temp;
		}

		chromosomalDistancesList = new ArrayList<Integer>(Arrays.asList(chromosomalDistancesArray));

		DistanceScore(nmerSize, resolution,type);

	}

	/*private void DistanceScore(int nmerSize, int resolution,String type) {
		//input
		String filepath = "output/" + type + "/" +Integer.toString(resolution) + "/" + Integer.toString(nmerSize) + "mersorted.txt";
		
		//output
		String outputFilepath = "output/" +type+"/"+ Integer.toString(resolution) + "/" + Integer.toString(nmerSize)
				+ "merscored.txt";

		
		int count = 0;
		try {
			BufferedReader bR = new BufferedReader(new FileReader(filepath));
			String str;
			String currentDistance = "";
			ArrayList<Integer> barcodeSizes = new ArrayList<Integer>();
			ArrayList<String> stringList = new ArrayList<String>();
			
			while ((str = bR.readLine()) != null) {
				count += 1;
				if (count % 100000 == 0) {
					System.out.println(count);
				}
				String[] splitLine = str.split(",");
				
				//keep looping through the sorted file until all k-mers with a given distance
				//are found and then score
				if (!currentDistance.equals(splitLine[2])) {
					writeScores(currentDistance, stringList, outputFilepath, barcodeSizes, resolution);
					currentDistance = splitLine[2];
					barcodeSizes = new ArrayList<Integer>();
					stringList = new ArrayList<String>();
				}
				barcodeSizes.add(Integer.parseInt(splitLine[1]));
				stringList.add(str);
			}
			bR.close();
		}

		catch (IOException e) {
			System.out.println("IO Exception! In getting score file");
			System.out.println(filepath);
		}

	}*/
	
	
	private void DistanceScore(int nmerSize, int resolution,String type) {
		//input
		String filepath = "output/" + type + "/" +Integer.toString(resolution) + "/" + Integer.toString(nmerSize) + "mersorted.txt";
		
		//output
		String outputFilepath = "output/" +type+"/"+ Integer.toString(resolution) + "/" + Integer.toString(nmerSize)
				+ "merscored.txt";

		
		int count = 0;
		try {
			BufferedReader bR = new BufferedReader(new FileReader(filepath));
			String str;
			String currentDistance = "";
			ArrayList<Integer> barcodeSizes = new ArrayList<Integer>();
			ArrayList<String> stringList = new ArrayList<String>();
			
			while ((str = bR.readLine()) != null) {
				count += 1;
				if (count % 100000 == 0) {
					System.out.println(count);
				}
				String[] splitLine = str.split(",");
				
				//String distance=splitLine[2];
				
				//keep looping through the sorted file until all k-mers with a given distance
				//are found and then score
				if (!currentDistance.equals(splitLine[2])) {
					writeScores(currentDistance, stringList, outputFilepath, barcodeSizes, resolution);
					currentDistance = splitLine[2];
					barcodeSizes = new ArrayList<Integer>();
					stringList = new ArrayList<String>();
				}
				barcodeSizes.add(Integer.parseInt(splitLine[1]));
				stringList.add(str);
			}
			writeScores(currentDistance, stringList, outputFilepath, barcodeSizes, resolution);
			bR.close();
		}

		catch (IOException e) {
			System.out.println("IO Exception! In getting score file");
			System.out.println(filepath);
		}

	}

	

	private void writeScores(String distance, ArrayList<String> strs, String outputFilepath, ArrayList<Integer> barcodeSizes, int resolution) throws IOException {
		FileWriter fw = new FileWriter(outputFilepath, true);

		BufferedWriter bw = new BufferedWriter(fw);
		PrintWriter out = new PrintWriter(bw);
		Collections.sort(barcodeSizes);
		
		//loop through all lines for a distance
		for (String str : strs) {
			String[] splitLine = str.split(",");
			int barcodeNumber = Integer.parseInt(splitLine[1]);
			long numberOfPossibleLocations = getPossibleLocations(distance, resolution);

			double averageBarcodeSize = getAverage(barcodeSizes);
			
			//compute proportion of k-mers for this distance greater than or equal 
			int totalGreaterThanOrEqual = barcodeSizes.size() - barcodeSizes.indexOf(barcodeNumber);

			double percentageRank = (double) totalGreaterThanOrEqual / numberOfPossibleLocations;

			out.println(String.join(",", splitLine) + "," + Double.toString(percentageRank) + ","
					+ Double.toString(barcodeNumber / averageBarcodeSize));
		}
		out.close();

	}


	private double getAverage(ArrayList<Integer> allBarcodeSizes) {
		int sum = 0;
		for (Integer size : allBarcodeSizes) {
			sum += size;
		}
		return sum / (double) allBarcodeSizes.size();
	}
	
	//get largest distances between all parts on the same chromsoome
	private ArrayList<Integer> parseMaxDistances(String distances) {
		ArrayList<Integer> returnDistances = new ArrayList<Integer>();
		String[] split = distances.split(":");

		for (int i = 0; i < split.length; i++) {
			if (split[i].indexOf("-") == -1) {
				returnDistances.add(Integer.parseInt(split[i]));
			} else {
				returnDistances.add(Integer.parseInt(split[i].substring(split[i].lastIndexOf("-") + 1)));
			}
		}
		return returnDistances;
	}

	private long getPossibleLocations(String distances, int resolution) {
		// get the maximum distance between points on the same chromosome with the kmers
		ArrayList<Integer> maxDistances = parseMaxDistances(distances);
		// nasty hack to have a long that can be passed by reference
		Long[] allLocationsSum = { 0L };
		// get all permutations of chromosomes of the same number as the number of
		// chromosomes in the kmer
		Generator.combination(chromosomalDistancesList).simple(maxDistances.size()).stream()
				.forEach(combination -> Generator.permutation(combination).simple().stream()
						.forEach(x -> appendDistances(x, maxDistances, allLocationsSum, resolution)));

		return allLocationsSum[0];
	}

	// for all combinations of chromosomes 
	private void appendDistances(List<Integer> chromosomeLength, ArrayList<Integer> maxDistances,
			Long[] allLocationsSum, int resolution) {
		long allCombinationsForChromosomePermutation = 1;
		for (int i = 0; i < maxDistances.size(); i++) {
			int allCombinationsForThisChromosome = (int) (Math.round(chromosomeLength.get(i) / (double) resolution)
					- maxDistances.get(i));
			// if region span too wide for this chromosome the permutation is invalid
			// and thus contributes no additional links
			if (allCombinationsForThisChromosome <= 0) {
				allCombinationsForChromosomePermutation = 0;
				break;
			}
			allCombinationsForChromosomePermutation *= allCombinationsForThisChromosome;
		}
		allLocationsSum[0] += allCombinationsForChromosomePermutation;

	}

	public static void main(String[] args) {
		System.out.println("Start distance scoring");
		for (int i = 0; i < args.length; i++) {
			System.out.println(args[i]);
		}
		int resolution = Integer.parseInt(args[0]);
		String type = args[1];
		
		for (int i = 2; i < args.length; i++) {

		new DistanceScoring(Integer.parseInt(args[i]), resolution, type);

		}
	}

}
