package guttmanlab.core.integrationTests;

import java.io.File;
import java.util.Random;

import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class OverlapperCountsIT {
    /*
     * Test code base counts against direct results from sam queries. Generates a random interval on chr19 and compares results
     * Compares results of BAMSingleReadCollection.numOverlappers() to SAMReader.queryOverlapping() or SAMReader.queryContained()
     */
    public static void main(String args[]) throws Exception{
        File bamfile = new File("/Users/cburghard/Downloads/chr19.clean.sorted.bam");
        BAMSingleReadCollection reads = new BAMSingleReadCollection(bamfile);

        int num_trials = 1000;
        Random rng = new Random();
        int wrongCount = 0;
        int c = 0;
        for(int i=0; i < num_trials;i++)
        {

            String refName = "chr19";
            int win_size = rng.nextInt(1000000);
            int start_pos = rng.nextInt(61342429-win_size); //chr19 size
            int end_pos = start_pos + rng.nextInt(win_size);

            int samOver = samOverlapping(refName,start_pos,end_pos);
            int samCont = samContained(refName,start_pos,end_pos);
            int over = guttmanOverlapping(reads,refName,start_pos,end_pos);
            int cont = guttmanContained(reads,refName,start_pos,end_pos);

            if(samOver != 0)
                c++;
            if(samCont != cont)
            {
                wrongCount++;
                System.out.println(refName + ":" + start_pos + "-" + end_pos);
                System.out.println("\tOverlap: " + samOver + "\t" + over);
                System.out.println("\tContained: " + samCont + "\t" + cont + "\n");
            }
        }
        System.out.println(wrongCount+"/"+c);
    }


    private static int samOverlapping(String refName,int start_pos,int end_pos)
    {
        File bamfile = new File("/Users/cburghard/Downloads/chr19.clean.sorted.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(bamfile);
        SAMRecordIterator currentIterator = reader.queryOverlapping(refName
                ,start_pos
                ,end_pos);

        int count = 0;
        while (currentIterator.hasNext())
        {
            currentIterator.next();
            count++;
        }

        return count;
    }

    private static int samContained(String refName,int start_pos,int end_pos)
    {
        File bamfile = new File("/Users/cburghard/Downloads/chr19.clean.sorted.bam");
        final SamReader reader = SamReaderFactory.makeDefault().open(bamfile);
        SAMRecordIterator currentIterator = reader.queryContained(refName
                ,start_pos
                ,end_pos);

        int count = 0;
        while (currentIterator.hasNext())
        {
            currentIterator.next();
            count++;
        }

        return count;
    }

    private static int guttmanOverlapping(BAMSingleReadCollection reads, String refName, int start_pos,int end_pos) throws Exception
    {
        return reads.numOverlappers(new SingleInterval(refName,start_pos,end_pos,Strand.BOTH), false);
    }

    private static int guttmanContained(BAMSingleReadCollection reads, String refName, int start_pos,int end_pos) throws Exception
    {
        return reads.numOverlappers(new SingleInterval(refName,start_pos,end_pos,Strand.BOTH), true);
    }
}
