package guttmanlab.core.integrationTests;

import guttmanlab.core.annotation.DerivedAnnotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SAMFragment;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.annotationcollection.AnnotationCollection;
import guttmanlab.core.annotationcollection.BAMSingleReadCollection;
import guttmanlab.core.coordinatespace.CoordinateSpace;
import net.sf.samtools.SAMFileHeader;

import java.io.File;
import java.io.IOException;

/**
 * Created by cburghard on 9/29/15.
 */
public class ConvertCoordinateSpaceIT {
    public static void main(String args[]) throws IOException {
        String pwd = args[0];
        BAMSingleReadCollection bam = new BAMSingleReadCollection(new File(pwd + args[1]));
        SAMFileHeader fhead = bam.getFileHeader();
        CoordinateSpace refSpace = new CoordinateSpace(fhead);
        String fname = pwd + "top50RefSeq.bed";
        BEDFileIO io = new BEDFileIO(pwd + "refspace.txt");
        AnnotationCollection<Gene> features = io.loadFromFile(fname);
        CoordinateSpace featureSpace = new CoordinateSpace(fname);

        AnnotationCollection<DerivedAnnotation<SAMFragment>> result = features.convertCoordinates(bam, refSpace, true);
        result.writeToBAM(pwd + args[2]);
    }
}
