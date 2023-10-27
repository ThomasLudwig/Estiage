package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Phase input data<br/>
 * 2 columns per marker (1 per chromosome)<br/>
 * 1 line per sample <br/>
 * 1 line header "m1 m1 m2 m2 .... mN mN
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-05-25
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Unphased {
  final int target;
  private final String EMPTY = "0";
  final String[] markerNames;
  final String[] sampleNames;
  final String[][][] inputData;
  final String[][] phasedData;
  String targetMarker = "TARGET";
  final int S;
  final int C;
  final String position;

  public Unphased(String filename, int target, String position) throws IOException, EstiageFormatException {
    this.target = target;
    this.position = position;
    UniversalReader in = new UniversalReader(filename);
    ArrayList<String[]> lines = new ArrayList<>();
    String line;
    while((line = in.readLine()) != null){
      lines.add(line.split("\t"));
    }
    in.close();
    String[] tmp = lines.get(0);
    String[] header = Arrays.copyOfRange(tmp, 1, tmp.length);
    S = lines.size() - 1;
    final int N = header.length;
    C = N  / 2;
    sampleNames = new String[S];

    //check parity
    if( 2*C != N )
      throw new EstiageFormatException("Expecting an even number of columns (2 per marker");

    //create markerNames
    markerNames = new String[C];
    //populate markerNames
    for(int i = 0 ; i < C; i++){
      markerNames[i] = header[i*2];
      if(!header[i*2 + 1].equals(markerNames[i]))
        throw new EstiageFormatException("Marker name mismatch for column ["+i+"]");
    }
    //check length of the rest of the line
    //create input
    inputData = new String[S][C][2];
    for(int s = 0; s < S; s++){
      String[] l = lines.get(s+1);
      sampleNames[s] = l[0];
      for(int i = 0 ; i < C; i++){
        inputData[s][i][0] = l[i*2 + 1];
        inputData[s][i][1] = l[i*2 + 2];
      }
    }
    //create output
    phasedData = new String[S][C];
    //populate input
  }

  public void phase() throws EstiageFormatException {
    System.err.println("Phasing for Marker ["+markerNames[target - 1]+"]");
    boolean[] keepLeft = new boolean[S];

    System.err.println("Left");
    //left
    for(int s = 0; s < S; s++)
      keepLeft[s] = true;
    rank(keepLeft, target-1);
    boolean[] keepRight = keepLeft.clone();
    for(int c = target - 2; c >= 0; c--) {
      System.err.println("Column "+c);
      for(int s = 0; s < S; s++)
        if(!keepLeft[s])
          phasedData[s][c] = "";
      if(rank(keepLeft, c)) {
        System.err.println("Ex aequo on column " + c + ", you should not go left from this point");
        break;
      }
    }

    //right
    for(int c = target; c < C; c++) {
      System.err.println("Column "+c);
      for (int s = 0; s < S; s++)
        if (!keepRight[s])
          phasedData[s][c] = "";
      if(rank(keepRight, c)) {
        System.err.println("Ex aequo on column " + c + ", you should not go right from this point");
        break;
      }
    }
  }

  private String getAllele(String[] ranks, String[] geno){
    for(String rank : ranks)
      if(rank.equals(geno[0]) || rank.equals(geno[1]))
        return rank;
    return EMPTY;
  }

  private boolean rank(boolean[] keep, int col) throws EstiageFormatException {
    HashMap<String, Integer> count = new HashMap<>();
    int dropped = 0;
    int empty = 0;
    for (int s = 0; s < S; s++)
      if (keep[s]) {
        if(inputData[s][col][0].isEmpty() || inputData[s][col][0].equals(EMPTY) || inputData[s][col][1].isEmpty() || inputData[s][col][1].equals(EMPTY)) {
          empty++;
          if(inputData[s][col][0].equals(EMPTY))
            inputData[s][col][0] = "";
          if(inputData[s][col][1].equals(EMPTY))
            inputData[s][col][1] = "";
        }
        else {
          Integer v1 = count.get(inputData[s][col][0]);
          if (v1 == null)
            v1 = 0;
          v1++;
          count.put(inputData[s][col][0], v1);
          if(!inputData[s][col][0].equals(inputData[s][col][1])){
            Integer v2 = count.get(inputData[s][col][1]);
            if (v2 == null)
              v2 = 0;
            v2++;
            count.put(inputData[s][col][1], v2);
          }
        }
      } else
        dropped++;

    System.err.println("dropped / empty : "+dropped+" / "+empty);
    int size = count.size();
    String[] ranks = new String[size];
    int[] vals = new int[size];

    for(int idx = 0; idx < size; idx++) {
      //Get allele with max value
      int max = 0;
      String allele = null;
      for (String key : count.keySet()) {
        if (count.get(key) > max) {
          max = count.get(key);
          allele = key;
        }
      }
      vals[idx] = max;
      ranks[idx] = allele;
      count.remove(allele);
    }

    System.err.println("For Marker ("+(col+1)+") ["+markerNames[col]+"]");
    for(int i = 0; i < size; i++)
      System.err.println("("+vals[i]+") -> "+ranks[i]);

    for(int s = 0; s < S; s++)
      if(keep[s]){
        //Apply allele to the previous position
        String allele = getAllele(ranks, inputData[s][col]);
        phasedData[s][col] = allele;
        if(!allele.equals(ranks[0])) {
          keep[s] = false;
        }
      }

    if(targetMarker.equals("TARGET"))
      targetMarker = ranks[0];
    return (size < 2 || vals[1] == vals[0]);
  }

  public void export(String outname) throws IOException {
    ArrayList<String> top = new ArrayList<>();
    int left = 0;
    int right = 0;
    top.add(".");

    ArrayList<String> line = new ArrayList<>();
    line.add("Samples\\Markers");
    for(int i = 0 ; i < markerNames.length; i++) {
      if(phasedData[0][i] != null) {
        line.add(markerNames[i]);
        if(i < target-1)
          left++;
        if(i >= target)
          right++;
      }
    }
    String pos = "Positions";
    for(int i = left; i > 0; i--) {
      top.add("Left" + i);
      pos+="\t???:???";
    }
    top.add("Target");
    pos += "\t"+position;
    for(int i = 1; i <= right; i++) {
      top.add("Right" + i);
      pos+="\t???:???";
    }

    PrintWriter out = new PrintWriter(new FileWriter(outname));
    out.println(String.join("\t", top));
    out.println(String.join("\t", line));
    for(int s = 0 ; s < S; s++) {
      if(phasedData[s][target - 1].equals(targetMarker)) {
        ArrayList<String> sline = new ArrayList<>();
        sline.add(sampleNames[s]);
        for (String geno : phasedData[s])
          if (geno != null)
            sline.add(geno);
        out.println(String.join("\t", sline));
      }
    }

    out.println(pos);
    out.close();
  }
}
