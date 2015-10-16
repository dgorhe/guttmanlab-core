package guttmanlab.core.integrationTests;

import java.io.File;

import guttmanlab.core.annotationcollection.BAMPairedFragmentCollection;

public class PairedFragmentTests {

	public static void main(String[] args) {
		File bamfile = new File("/Users/cburghard/Downloads/chr19.clean.sorted.bam");
		BAMPairedFragmentCollection fragments = new BAMPairedFragmentCollection(bamfile);
		
		fragments.writeToBAM("/Users/cburghard/Downloads/chr19.fragments.bam");
	}

}
