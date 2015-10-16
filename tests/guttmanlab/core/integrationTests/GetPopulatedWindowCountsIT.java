package guttmanlab.core.integrationTests;

import guttmanlab.core.annotation.*;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AbstractAnnotationCollection;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMFragmentCollectionFactory;
import net.sf.samtools.util.CloseableIterator;

import java.io.*;
import java.util.ArrayList;

/**
 * Created by cburghard on 9/21/15.
 */
public class GetPopulatedWindowCountsIT {
    public static void main(String args[]) throws IOException{
        //1. test getWindows returns correct intervals over genes
        //10	NM_001163328	chr3	+	13471654	14182287	13776847	14181055	9
        //13471654,13776824,14039713,14081207,14107218,14115692,14129557,14143273,14181037,
        // 13471943,13777103,14039789,14081240,14107266,14115856,14129671,14143446,14182287

        BufferedReader br = new BufferedReader(new FileReader("/Users/cburghard/Downloads/test.bed"));
        String line;
        while ((line = br.readLine()) != null) {
            print_exons(line);
        }

        //get populated windows
        File bamFile = new File("/Users/cburghard/Documents/GuttmanLab/Ribosome19.bam");
        AbstractAnnotationCollection<? extends MappedFragment> bam= BAMFragmentCollectionFactory.createFromBam(bamFile, false);
        BEDFileIO io =  new BEDFileIO("/Users/cburghard/Documents/GuttmanLab/References/refspace.txt");
        AnnotationCollection<Gene> genes = io.loadFromFile("/Users/cburghard/Downloads/test.bed");

        CloseableIterator<Gene> g_iter = genes.sortedIterator();
        while( g_iter.hasNext())
        {
            Gene g = g_iter.next();
            AnnotationCollection<DerivedAnnotation<? extends Annotation>> windows = g.getWindows(100, 100);
            CloseableIterator<DerivedAnnotation<? extends Annotation>> w_iter = windows.sortedIterator();
            while(w_iter.hasNext())
            {
                Annotation w = w_iter.next();
                System.out.println(w.toBED());
            }
        }
    }

    //print out window coordinates
    //compare to correct values

    //Alternate, independent method for determining correct window intervals to test against Annotation.getWindows()
    private static void print_exons(String bedLine){
        String chr_name = bedLine.split("\t")[0];
        String[] sizes = bedLine.split("\t")[10].split(",");
        String[] starts = bedLine.split("\t")[11].split(",");
        ArrayList<Integer> exon_ends = new ArrayList<Integer>();
        ArrayList<Integer> exon_starts = new ArrayList<Integer>();

        for( int i = 0; i < starts.length; i++ )
        {
            exon_starts.add(i,Integer.parseInt(starts[i]));
            exon_ends.add(exon_starts.get(i) + Integer.parseInt(sizes[i]));
            //System.out.println(exon_starts.get(i) + "\t" + exon_ends.get(i) + "\t" + sizes[i]);
        }

        int win_size = 100;
        Annotation window = new BlockedAnnotation();
        int exon_end = 0;
        int exon_start = 0;
        int win_end = 0;

        int remaining = 0;
        for ( int i = 0; i < exon_starts.size(); i++ )
        {
            exon_start = exon_starts.get(i);
            exon_end = exon_ends.get(i);
            int win_start =  exon_start;

            while( (win_end = win_start + remaining) <= exon_end )
            {
                System.out.println(chr_name + "\t" +win_start+"\t"+win_end);
                win_start = win_end;
                remaining = win_size;
            }

            //left over length at end of exon
            remaining = win_end - exon_end;
            if(remaining % win_size != 0)
                System.out.println(chr_name + "\t" + win_start+"\t"+exon_end);

        }
    }
}