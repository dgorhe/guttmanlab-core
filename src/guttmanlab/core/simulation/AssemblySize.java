package guttmanlab.core.simulation;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 * A convenience class containing static maps of assembly sizes. Used to
 * construct <code>CoordinateSpace</code> objects.
 * <p>
 * This class does not handle unlocalized sequences, unplaced sequences, or
 * alternate loci.
 */
public final class AssemblySize {
    
    public static final Map<String, Integer> MM9;
    public static final Map<String, Integer> MM10;
    public static final Map<String, Integer> HG19;
    public static final Map<String, Integer> HG38;
    
    private AssemblySize() { }

    static {
        Map<String, Integer> mm9Map = new TreeMap<>();
        mm9Map.put("chr1",197195432);
        mm9Map.put("chr2",181748087);
        mm9Map.put("chr3",159599783);
        mm9Map.put("chr4",155630120);
        mm9Map.put("chr5",152537259);
        mm9Map.put("chr6",149517037);
        mm9Map.put("chr7",152524553);
        mm9Map.put("chr8",131738871);
        mm9Map.put("chr9",124076172);
        mm9Map.put("chr10",129993255);
        mm9Map.put("chr11",121843856);
        mm9Map.put("chr12",121257530);
        mm9Map.put("chr13",120284312);
        mm9Map.put("chr14",125194864);
        mm9Map.put("chr15",103494974);
        mm9Map.put("chr16",98319150);
        mm9Map.put("chr17",95272651);
        mm9Map.put("chr18",90772031);
        mm9Map.put("chr19",61342430);
        mm9Map.put("chrX",166650296);
        mm9Map.put("chrY",15902555);
        mm9Map.put("chrM",16299);
        mm9Map.put("chr13_random",400311);
        mm9Map.put("chr16_random",3994);
        mm9Map.put("chr17_random",628739);
        mm9Map.put("chr1_random",1231697);
        mm9Map.put("chr3_random",41899);
        mm9Map.put("chr4_random",160594);
        mm9Map.put("chr5_random",357350);
        mm9Map.put("chr7_random",362490);
        mm9Map.put("chr8_random",849593);
        mm9Map.put("chr9_random",449403);
        mm9Map.put("chrUn_random",5900358);
        mm9Map.put("chrX_random",1785075);
        mm9Map.put("chrY_random",58682461);
        
        MM9 = Collections.unmodifiableMap(mm9Map);
    }
    
    static {
        Map<String, Integer> mm10Map = new TreeMap<>();
        mm10Map.put("chr1", 195471971);
        mm10Map.put("chr2", 182113224);
        mm10Map.put("chr3", 160039680);
        mm10Map.put("chr4", 156508116);
        mm10Map.put("chr5", 151834684);
        mm10Map.put("chr6", 149736546);
        mm10Map.put("chr7", 145441459);
        mm10Map.put("chr8", 129401213);
        mm10Map.put("chr9", 124595110);
        mm10Map.put("chr10", 130694993);
        mm10Map.put("chr11", 122082543);
        mm10Map.put("chr12", 120129022);
        mm10Map.put("chr13", 120421639);
        mm10Map.put("chr14", 124902244);
        mm10Map.put("chr15", 104043685);
        mm10Map.put("chr16", 98207768);
        mm10Map.put("chr17", 94987271);
        mm10Map.put("chr18", 90702639);
        mm10Map.put("chr19", 61431566);
        mm10Map.put("chrX", 171031299);
        mm10Map.put("chrY", 91744698);
        //mm10Map.put("chrM", 16299);
        MM10 = Collections.unmodifiableMap(mm10Map);
    }
    
    static {
        Map<String, Integer> hg19Map = new TreeMap<>();
        
        hg19Map.put("chr1",249250621);
        hg19Map.put("chr2",243199373);
        hg19Map.put("chr3",198022430);
        hg19Map.put("chr4",191154276);
        hg19Map.put("chr5",180915260);
        hg19Map.put("chr6",171115067);
        hg19Map.put("chr7",159138663);
        hg19Map.put("chrX",155270560);
        hg19Map.put("chr8",146364022);
        hg19Map.put("chr9",141213431);
        hg19Map.put("chr10",135534747);
        hg19Map.put("chr11",135006516);
        hg19Map.put("chr12",133851895);
        hg19Map.put("chr13",115169878);
        hg19Map.put("chr14",107349540);
        hg19Map.put("chr15",102531392);
        hg19Map.put("chr16",90354753);
        hg19Map.put("chr17",81195210);
        hg19Map.put("chr18",78077248);
        hg19Map.put("chr20",63025520);
        hg19Map.put("chrY",59373566);
        hg19Map.put("chr19",59128983);
        hg19Map.put("chr22",51304566);
        hg19Map.put("chr21",48129895);
        hg19Map.put("chr6_ssto_hap7",4928567);
        hg19Map.put("chr6_mcf_hap5",4833398);
        hg19Map.put("chr6_cox_hap2",4795371);
        hg19Map.put("chr6_mann_hap4",4683263);
        hg19Map.put("chr6_apd_hap1",4622290);
        hg19Map.put("chr6_qbl_hap6",4611984);
        hg19Map.put("chr6_dbb_hap3",4610396);
        hg19Map.put("chr17_ctg5_hap1",1680828);
        hg19Map.put("chr4_ctg9_hap1",590426);
        hg19Map.put("chr1_gl000192_random",547496);
        hg19Map.put("chrUn_gl000225",211173);
        hg19Map.put("chr4_gl000194_random",191469);
        hg19Map.put("chr4_gl000193_random",189789);
        hg19Map.put("chr9_gl000200_random",187035);
        hg19Map.put("chrUn_gl000222",186861);
        hg19Map.put("chrUn_gl000212",186858);
        hg19Map.put("chr7_gl000195_random",182896);
        hg19Map.put("chrUn_gl000223",180455);
        hg19Map.put("chrUn_gl000224",179693);
        hg19Map.put("chrUn_gl000219",179198);
        hg19Map.put("chr17_gl000205_random",174588);
        hg19Map.put("chrUn_gl000215",172545);
        hg19Map.put("chrUn_gl000216",172294);
        hg19Map.put("chrUn_gl000217",172149);
        hg19Map.put("chr9_gl000199_random",169874);
        hg19Map.put("chrUn_gl000211",166566);
        hg19Map.put("chrUn_gl000213",164239);
        hg19Map.put("chrUn_gl000220",161802);
        hg19Map.put("chrUn_gl000218",161147);
        hg19Map.put("chr19_gl000209_random",159169);
        hg19Map.put("chrUn_gl000221",155397);
        hg19Map.put("chrUn_gl000214",137718);
        hg19Map.put("chrUn_gl000228",129120);
        hg19Map.put("chrUn_gl000227",128374);
        hg19Map.put("chr1_gl000191_random",106433);
        hg19Map.put("chr19_gl000208_random",92689);
        hg19Map.put("chr9_gl000198_random",90085);
        hg19Map.put("chr17_gl000204_random",81310);
        hg19Map.put("chrUn_gl000233",45941);
        hg19Map.put("chrUn_gl000237",45867);
        hg19Map.put("chrUn_gl000230",43691);
        hg19Map.put("chrUn_gl000242",43523);
        hg19Map.put("chrUn_gl000243",43341);
        hg19Map.put("chrUn_gl000241",42152);
        hg19Map.put("chrUn_gl000236",41934);
        hg19Map.put("chrUn_gl000240",41933);
        hg19Map.put("chr17_gl000206_random",41001);
        hg19Map.put("chrUn_gl000232",40652);
        hg19Map.put("chrUn_gl000234",40531);
        hg19Map.put("chr11_gl000202_random",40103);
        hg19Map.put("chrUn_gl000238",39939);
        hg19Map.put("chrUn_gl000244",39929);
        hg19Map.put("chrUn_gl000248",39786);
        hg19Map.put("chr8_gl000196_random",38914);
        hg19Map.put("chrUn_gl000249",38502);
        hg19Map.put("chrUn_gl000246",38154);
        hg19Map.put("chr17_gl000203_random",37498);
        hg19Map.put("chr8_gl000197_random",37175);
        hg19Map.put("chrUn_gl000245",36651);
        hg19Map.put("chrUn_gl000247",36422);
        hg19Map.put("chr9_gl000201_random",36148);
        hg19Map.put("chrUn_gl000235",34474);
        hg19Map.put("chrUn_gl000239",33824);
        hg19Map.put("chr21_gl000210_random",27682);
        hg19Map.put("chrUn_gl000231",27386);
        hg19Map.put("chrUn_gl000229",19913);
        hg19Map.put("chrM",16571);
        hg19Map.put("chrUn_gl000226",15008);
        hg19Map.put("chr18_gl000207_random",4262);
        
       
       
        
        HG19 = Collections.unmodifiableMap(hg19Map);
    }
    
    
    static {
        Map<String, Integer> hg38Map = new TreeMap<>();
        
        hg38Map.put("chr1",248956422);
        hg38Map.put("chr2",242193529);
        hg38Map.put("chr3",198295559);
        hg38Map.put("chr4",190214555);
        hg38Map.put("chr5",181538259);
        hg38Map.put("chr6",170805979);
        hg38Map.put("chr7",159345973);
        hg38Map.put("chrX",156040895);
        hg38Map.put("chr8",145138636);
        hg38Map.put("chr9",138394717);
        hg38Map.put("chr10",133797422);
        hg38Map.put("chr11",135086622);
        hg38Map.put("chr12",133275309);
        hg38Map.put("chr13",114364328);
        hg38Map.put("chr14",107043718);
        hg38Map.put("chr15",101991189);
        hg38Map.put("chr16",90338345);
        hg38Map.put("chr17",83257441);
        hg38Map.put("chr18",80373285);
        hg38Map.put("chr20",64444167);
        hg38Map.put("chrY",57227415);
        hg38Map.put("chr19",58617616);
        hg38Map.put("chr22",50818468);
        hg38Map.put("chr21",46709983);
        
        
    
        
       
       
        
        HG38 = Collections.unmodifiableMap(hg38Map);
    }
}