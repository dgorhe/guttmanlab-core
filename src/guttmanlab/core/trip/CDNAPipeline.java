package guttmanlab.core.trip;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import guttmanlab.core.annotation.io.BEDFileIO;

public class CDNAPipeline {

	public static void pipeline(String file, int distance) throws IOException {
		//Step 1: collapse barcodes --> print barcodes to UMIs
		List<String> lines=BEDFileIO.loadLines(file);
		Map<String, Map<String, String>> map=CollapseCDNA.collapseBarcodes(lines, distance);
		//List<String> files=write(map, save);
		
		System.err.println("barcode sizes "+map.keySet());
		
		
		//Step 2: parallelize for each barcode collapse by UMI
		for(String barcode: map.keySet()) {
			Collection<String> collapsedUMI=CollapseCDNA.collapseUMI(map.get(barcode),1);
			System.out.println(barcode+" "+map.get(barcode).size()+" "+collapsedUMI.size());
		}
		
		
		
	}

	
	
	
	public static void main(String[] args) throws IOException {
		pipeline(args[0], Integer.parseInt(args[1]));
	}
	
}
