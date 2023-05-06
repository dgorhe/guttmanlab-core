package guttmanlab.core.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

public class TestCoPilot {

    public static void main(String[] args) {
        System.out.println("Hello World!");


    }

    //a method that takes a collection of strings and collapses all strings that are within a certain edit distance of each other
    private void collapseStrings(Collection<String> input, int editDistance){
        input.size();
        input.iterator();


        //for each string in the input collection
        for(String s : input){
            //for each other string in the input collection
            for(String t : input){
                //if the edit distance between s and t is less than or equal to the edit distance argument
                if(editDistance(s, t) <= editDistance){
                    //remove t from the input collection
                    input.remove(t);
                }
            }
        }
    }

    //a method to compute the edit distance between two strings
    private int editDistance(String s, String t){
        //if the first string is empty, the edit distance is the length of the second string
        if(s.isEmpty()){
            return t.length();
        }
        //if the second string is empty, the edit distance is the length of the first string
        if(t.isEmpty()){
            return s.length();
        }
        //if the first character of the first string is the same as the first character of the second string
        if(s.charAt(0) == t.charAt(0)){
            //the edit distance is the edit distance between the rest of the first string and the rest of the second string
            return editDistance(s.substring(1), t.substring(1));
        }
        //otherwise the edit distance is the minimum of the edit distance between the rest of the first string and the second string,
        //the first string and the rest of the second string, and the rest of the first string and the rest of the second string
        return 1 + Math.min(editDistance(s.substring(1), t), Math.min(editDistance(s, t.substring(1)), editDistance(s.substring(1), t.substring(1))));


    }

    //a method to count the number of lines containing the same barcode string present in the first column of the file
    //import all relevent packages
    private void process(File file) throws IOException {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line;
        Map<String, Integer> barcodeCounts = new HashMap<String, Integer>();
        while((line = reader.readLine()) != null){
            String[] tokens = line.split("\t");
            String barcode = tokens[0];
            if(barcodeCounts.containsKey(barcode)){
                barcodeCounts.put(barcode, barcodeCounts.get(barcode) + 1);
            }else{
                barcodeCounts.put(barcode, 1);
            }
        }
        reader.close();
        for(String barcode : barcodeCounts.keySet()){
            System.out.println(barcode + "\t" + barcodeCounts.get(barcode));
        }
    }

}
