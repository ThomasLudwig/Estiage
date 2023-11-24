package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * Phase input data<br/>
 * 2 columns per marker (1 per chromosome)<br/>
 * 1 header line "m1a m1b m2a m2b ... T1a T2b ... mNa mNb" (where T is the Target Marker)
 * 1 line per sample <br/>
 * no additional lines
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-05-25
 * Checked for release on 2023-11-24
 * Unit Test defined on   2023-11-24
 */
public class Unphased {
  public static final String LEFT = "left";
  public static final String RIGHT = "Right";
  public static final String EMPTY = "0";
  private static final String UNDEFINED_ALLELE = "UNDEFINED_ALLELE";

  final int target;
  final String[] markerNames;
  final String[] sampleNames;
  final String[][][] inputData;

  final int S;
  final int C;
  final String position;

  /**
   * Loads an Unphased input file and create a phased version
   * @param filename the name of the input file
   * @param target the index of the column containing the target marker
   * @param position the position of the target marker (chr:pos)
   * @throws IOException if the input file can't be read
   * @throws EstiageFormatException if the input file can't be parsed
   */
  public Unphased(String filename, int target, String position) throws IOException, EstiageFormatException {
    this.target = target;
    this.position = position;

    //Load all the lines
    UniversalReader in = new UniversalReader(filename);
    ArrayList<String[]> lines = new ArrayList<>();
    String line;
    while((line = in.readLine()) != null){
      lines.add(line.split("\t", -1));
    }
    in.close();

    //Split first line and ignore first column
    String[] tmp = lines.get(0);
    Message.debug("There are "+tmp.length+" columns in the input file");
    String[] header = Arrays.copyOfRange(tmp, 1, tmp.length);

    //Get the number of samples (number of lines minus header)
    S = lines.size() - 1;
    //Get the number of position/markers (number of non header column / 2)
    final int N = header.length;
    C = N  / 2;
    Message.debug("There are "+C+" markers and "+S+" samples");
    sampleNames = new String[S];

    //check parity
    if( 2*C != N )
      throw new EstiageFormatException("Expecting an even number of columns (2 per marker");

    //create markerNames
    markerNames = new String[C];
    //populate markerNames
    for(int i = 0 ; i < C; i++){
      markerNames[i] = header[i*2];
      if(!header[i*2].equals(markerNames[i]))
        throw new EstiageFormatException("Marker name mismatch for column ["+i+"]");
    }

    //create input
    inputData = new String[S][C][2];
    for(int s = 0; s < S; s++){
      String[] l = lines.get(s+1);
      //check length of the line
      if(l.length != N+1)
        throw new EstiageFormatException("For line ["+(s+1)+"] the number of columns found ["+l.length+"] differs from the number of columns expected ["+(N+1)+"]");
      //populate input
      sampleNames[s] = l[0];
      for(int i = 0 ; i < C; i++){
        inputData[s][i][0] = l[i*2 + 1];
        inputData[s][i][1] = l[i*2 + 2];
      }
    }
  }

  /**
   * The phased version of the data
   */
  static class Phased {
    private final String[] sampleNames;
    private final String[] markerNames;
    private final String[][][] input;
    private final String[][] data;
    private final int target;
    private final String position;
    private String targetMarkerAllele;
    private final int S;
    private final int C;
    private final boolean stopOnExAequo;
    private final boolean ignoreMissing;

    /**
     * Build phased data from unphased ones
     * @param unphased the source unphased data
     */
    public Phased(Unphased unphased, boolean stopOnExAequo, boolean ignoreMissing){
      this.sampleNames = unphased.sampleNames.clone();
      this.markerNames = unphased.markerNames.clone();
      this.input = unphased.inputData.clone();
      this.position = unphased.position;
      this.targetMarkerAllele = UNDEFINED_ALLELE;
      this.target = unphased.target;
      S = input.length;
      C = input[0].length;
      data = new String[S][C];
      this.stopOnExAequo = stopOnExAequo;
      this.ignoreMissing = ignoreMissing;
      phase();
    }

    /**
     * Creates phased data out of the unphased ones
     */
    public void phase() {
      Message.info("Phasing for Marker ["+markerNames[target - 1]+"]");
      //init keep is regard to target
      //kept samples on the left
      boolean[] keepLeft = new boolean[S];
      for(int s = 0; s < S; s++)
        keepLeft[s] = true;
      rank(keepLeft, target-1);
      //kept samples on the right
      boolean[] keepRight = keepLeft.clone();

      //left
      Message.info("Left");
      boolean stop = false;
      int last = target - 1;
      for(int c = target - 2; c >= 0 && !stop; c--) {
        stop = phasePosition(c, keepLeft, LEFT);
        last = c;
      }
      trackBackExAequo(last, target - 1);
      //right
      Message.info("Right");
      stop = false;
      for(int c = target; c < C && !stop; c++) {
        stop = phasePosition(c, keepRight, RIGHT);
        last = c;
      }
      trackBackExAequo(last, target - 1);
    }

    /**
     *
     * @param c the current marker index (column)
     * @param keep the list of individual to keep(true)/remove(false)
     * @param side the current side's name
     * @return false if the algorithm needs to stop for this side
     */
    private boolean phasePosition(int c, boolean[] keep, String side){
      Message.info("Column "+c);
      for (int s = 0; s < S; s++)
        if (!keep[s])
          data[s][c] = "";
      boolean stop = rank(keep, c);
      if(stop)
        Message.info(c + " is the last column on the "+side);
      return stop;
    }

    /**
     * Replace default "empty" value with an empty string and check if at list one chromosome is empty
     * @param val the genotype to consider
     * @return true is at least one allele is empty
     */
    private static boolean isEmpty(String[] val) {
      for(int i = 0 ; i < 2 ; i++)
        if(EMPTY.equals(val[i]))
          val[i] = "";
      if(val[0].isEmpty() != val[1].isEmpty())
        Message.warning("Only one genotype is empty and not the other");

      return val[0].isEmpty() && val[1].isEmpty() ;
    }

    /**
     *  1. Remove samples that do not carry the top allele
     *  2. affect top allele to phased data if the sample carries it
     *  3. Reset the target market value if needed
     *  4. check if unique value or ex aequo
     * @param keep the samples to keep (array of bool. true at ith index means the ith individual is kept)
     * @param col the column (index of marker) to consider
     * @return true if algorithm needs to stop for this side
     */
    private boolean rank(boolean[] keep, int col) {
      HashMap<String, Integer> counts = new HashMap<>();
      int dropped = 0;
      int empty = 0;

      //for each sample, count each genotype (only once for homozygous)
      for (int s = 0; s < input.length; s++) {
        //if the sample is still kept
        if (keep[s]) {
          //if there is an empty value on either chromosome per the current sample/marker
          if (isEmpty(input[s][col])) {
            empty++;
          } else {
            //increment count for distinct alleles
            HashSet<String> distinctAlleles = new HashSet<>();
            for(String allele: input[s][col])
              if(!allele.isEmpty())
                distinctAlleles.add(allele);
            for(String allele : distinctAlleles)
              counts.merge(allele, 1, Integer::sum);
          }
        } else
          dropped++;
      }

      Message.info("dropped / empty : "+dropped+" / "+empty);
      Ranking<String> ranking = new Ranking<>(counts);
      String ancestral = ranking.getElement(0);

      Message.debug("For Marker (" + (col + 1) + ") [" + markerNames[col] + "]");
      for(int i = 0; i < ranking.size(); i++)
        Message.debug("("+ranking.getOccurrence(i)+") -> "+ranking.getElement(i));

      boolean stop = true;
      // for each kept sample
      for(int s = 0; s < S; s++) {
        if (keep[s]) {
          //Apply allele to the previous position
          Genotype geno = new Genotype(ranking, input[s][col], stopOnExAequo, ignoreMissing);
          data[s][col] = geno.toString();
          keep[s] = geno.isKept();
          if(!keep[s])
            Message.debug("Remove " + sampleNames[s] + " [" + geno + "]");
        }
        if(keep[s])
          stop = false;
      }

      if(targetMarkerAllele.equals(UNDEFINED_ALLELE))
        targetMarkerAllele = ancestral;

      if(stop)
        return true;

      if(ignoreMissing)
        return ranking.size() < 2;
      return (ranking.size() < 2 || ranking.hasTopExAequo());
    }

    private void trackBackExAequo(int from, int to){
      //TODO track back ex aequo at previous rank, if a current marker clear a previous ex aequo
      //For the moment it is to much work for something that could be done quickly by a human.
      //It isn't even clear if
      // - it is always possible
      // - we should start from the extremity or the target

      //step one
      //  select all the sample with the longest haplotype
      //  check if at any position they carry an allele that is ex aequo, and lift this ex aequo


      int p = 1;
      if(to < from)
        p = -1;
      for(int c = from; c != to; c+=p){

      }
    }

    /**
     * Exports the phased version of the input file
     * @param filename the name of the output file
     * @throws IOException if the file can't be written
     */
    public void export(String filename) throws IOException {
      ArrayList<String> header1 = new ArrayList<>();//. LeftN...Left2 Left1 Target Right1 Right2...RightN
      ArrayList<String> header2 = new ArrayList<>();//Samples\Markers D1S121...D1S417 D1S478 D1S587 D1S237...D1S334
      header1.add(".");
      header2.add("Samples\\Markers");

      int left = 0;
      int right = 0;

      //Fill headers and Position
      for(int i = 0 ; i < markerNames.length; i++) {
        if(data[0][i] != null) {
          header2.add(markerNames[i]);
          if(i < target-1)
            left++;
          if(i >= target)
            right++;
        }
      }
      StringBuilder pos = new StringBuilder("Positions");
      for(int i = left; i > 0; i--) {
        header1.add("Left" + i);
        pos.append("\t???:???");
      }
      header1.add("Target");
      pos.append("\t").append(position);
      for(int i = 1; i <= right; i++) {
        header1.add("Right" + i);
        pos.append("\t???:???");
      }

      //Print the headers, data and position
      PrintWriter out = new PrintWriter(new FileWriter(filename));
      out.println(String.join("\t", header1));
      out.println(String.join("\t", header2));
      for(int s = 0 ; s < data.length; s++) {
        if(data[s][target - 1].equals(targetMarkerAllele)) { //ignore sample without target allele
          ArrayList<String> sLine = new ArrayList<>();
          sLine.add(sampleNames[s]);
          for (String geno : data[s])
            if (geno != null)
              sLine.add(geno);
          out.println(String.join("\t", sLine));
        }
      }
      out.print(pos);
      out.flush();
      out.close();
    }
  }

  private static class Genotype {
    private final String allele1;
    private final String allele2;
    private final boolean keep;
    private final boolean stopOnExAequo;
    private final boolean ignoreMissing;

    public Genotype(Ranking<String> ranking, String[] geno, boolean stopOnExAequo, boolean ignoreMissing){
      this.stopOnExAequo = stopOnExAequo;
      this.ignoreMissing = ignoreMissing;
      String[] a = {EMPTY, EMPTY};
      int[] o = {-1,-1};

      for(int i = 0 ; i < ranking.size(); i++) {
        for(int g = 0 ; g < 2; g++)
          if(geno[g].equals(ranking.getElement(i))){
            a[g] = geno[g];
            o[g] = ranking.getOccurrence(i);
          }
      }

      if(o[0] > o[1]){
        this.allele1 = a[0];
        this.allele2 = this.allele1;
      } else if(o[1] > o[0]) {
        this.allele1 = a[1];
        this.allele2 = this.allele1;
      } else {
        this.allele1 = min(a[0], a[1]);
        this.allele2 = max(a[0], a[1]);
      }
      this.keep = this.computeKeep(ranking);
    }

    private static String min(String s1, String s2){
      if(s1.compareTo(s2) <= 0)
        return s1;
      return s2;
    }

    private static String max(String s1, String s2){
      if(s1.compareTo(s2) <= 0)
        return s2;
      return s1;
    }

    public boolean isEmpty(){
      return EMPTY.equals(allele1) && EMPTY.equals(allele2);
    }

    public boolean isKept(){
      return this.keep;
    }

    private boolean computeKeep(Ranking<String> ranking){
      if (stopOnExAequo){
        if (ranking.hasTopExAequo())
          return false;
      }
      if (ignoreMissing){
        if(isEmpty())
          return true;
      }
      return ranking.getTopElements().contains(this.allele1);
    }

    @Override
    public String toString() {
      if(this.allele1.equals(allele2))
        return allele1;
      return allele1+"/"+allele2;
    }
  }
}
