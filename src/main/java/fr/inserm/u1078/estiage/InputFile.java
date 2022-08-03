package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Estiage Original Input File Format
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2021-03-23
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class InputFile {
  private final int nbSamples;
  private final double mutationRate;
  private final int mutationModel;
  /*
    private final int nbLeft;
    private final int nbRight;
    private final double[] leftFractions;
    private final double[] rightFractions;
    private final double[] leftFrequencies;
    private final double[] rightFrequencies;
    private final String[] leftAncestral;
    private final String[] rightAncestral;
    private final String[][] records;
  */

  private final Side left;
  private final Side right;

  public InputFile(TSVFile raw, int mutationModel, double mutationRate) {
    this.nbSamples = raw.getSamples().length;
    this.mutationRate = mutationRate;
    this.mutationModel = mutationModel;
    /*this.nbLeft = raw.getLeftMarkers().length;
    this.nbRight = raw.getRightMarkers().length;
    this.leftFractions = new double[nbLeft];
    this.leftFrequencies = new double[nbLeft];
    this.leftAncestral = new String[nbLeft];
    this.rightFractions = new double[nbRight];
    this.rightFrequencies = new double[nbRight];
    this.rightAncestral = new String[nbRight];

    this.records = new String[nbSamples][4];
    for (int i = 0; i < nbSamples; i++)
      for (int j = 0; j < 4; j++)
        this.records[i][j] = "-1";*/

    left = load(/*LEFT, */raw.getLeftMarkers(), raw.getSamples());//, leftFractions, leftFrequencies, leftAncestral);
    right = load(/*RIGHT, */raw.getRightMarkers(), raw.getSamples());//, rightFractions, rightFrequencies, rightAncestral);
  }

  /**
   * Constructor. Builds an InputFile object from a serialized input file
   * @param preinputFile
   * @throws EstiageFormatException
   * @throws IOException
   */
  public InputFile(String preinputFile) throws EstiageFormatException, IOException {
    UniversalReader in = new UniversalReader(preinputFile);
    String[] f;
    int line = 0;
    try {
      //header
      line++;
      f = in.readLine().split(S, -1);
      this.nbSamples = new Integer(f[0]);
      int nbLeft = new Integer(f[1]);
      int nbRight = new Integer(f[2]);

      double[] fractionsLeft = new double[nbLeft];
      String[] frequenciesLeft = new String[nbLeft];
      String[] ancestralLeft = new String[nbLeft];

      double[] fractionsRight = new double[nbRight];
      String[] frequenciesRight = new String[nbRight];
      String[] ancestralRight = new String[nbRight];

      int[] positionsLeft = new int[nbSamples];
      String[] endLeft = new String[nbSamples];
      int[] positionsRight = new int[nbSamples];
      String[] endRight = new String[nbSamples];

      //fracleft
      line++;
      f = in.readLine().split(S, -1);
      for(int i = 0 ; i < nbLeft; i++)
        fractionsLeft[i] = new Double(f[i]);

      //freqleft
      line++;
      f = in.readLine().split(S, -1);
      for(int i = 0 ; i < nbLeft; i++)
        frequenciesLeft[i] = f[i];

      //fracright
      line++;
      f = in.readLine().split(S, -1);
      for(int i = 0 ; i < nbRight; i++)
        fractionsRight[i] = new Double(f[i]);

      //freqright
      line++;
      f = in.readLine().split(S, -1);
      for(int i = 0 ; i < nbRight; i++)
        frequenciesRight[i] = f[i];

      //mutation
      line++;
      f = in.readLine().split(S, -1);
      this.mutationRate = new Double(f[0]);
      this.mutationModel = new Integer(f[1]);

      //ancestralleft
      line++;
      f = in.readLine().split(S, -1);
      for(int i = 0 ; i < nbLeft; i++)
        ancestralLeft[i] = f[i];

      //ancestralright
      line++;
      f = in.readLine().split(S, -1);
      for(int i = 0 ; i < nbRight; i++)
        ancestralRight[i] = f[i];

      // samples
      for(int i = 0 ; i < nbSamples; i++) {
        line++;
        f = in.readLine().split(S, -1);
        positionsLeft[i] = new Integer(f[0]);
        positionsRight[i] = new Integer(f[1]);
        endLeft[i] = f[2];
        endRight[i] = f[3];
      }

      this.left = new Side(nbLeft, fractionsLeft, frequenciesLeft, ancestralLeft, positionsLeft, endLeft);
      this.right = new Side(nbRight, fractionsRight, frequenciesRight, ancestralRight, positionsRight, endRight);
    } catch(NullPointerException | NumberFormatException | ArrayIndexOutOfBoundsException e) {
      throw new EstiageFormatException("Unable to parse input file ["+preinputFile+"] on line ["+line+"]", e);
    }
    in.close();
  }

  /**
   * Converts a preinput file to an input file
   */
  public void fromPreinput2Input(){
    left.fromPreinput2Input();
    right.fromPreinput2Input();
  }

  public static boolean noValue(String s){
    if("-1".equals(s))
      return true;
    return s.isEmpty();
  }

  /**
   * Creates a Side object from Markers and Samples
   * @param markers
   * @param samples
   * @return
   */
  private Side load(Marker[] markers, String[] samples){
    int nb = markers.length;
    boolean ambiguous = !noValue(markers[nb-1].getAncestral());

    if(ambiguous)
      nb++;

    double[] fractions = new double[nb];
    String[] frequencies = new String[nb];
    String[] ancestrals = new String[nb];

    for (int i = 0; i < markers.length; i++) {
      Marker m = markers[i];
      String ancestral = m.getAncestral();
      fractions[i] = m.getRecombinationFraction();
      frequencies[i] = ""+m.getFrequencies();
      ancestrals[i] = ancestral;
    }

    if(ambiguous){//If the last marker has a clear ancestral allele, add another marker with a fraction of .5 after this one
      fractions[nb - 1] = .5;
      frequencies[nb - 1] = "";
      ancestrals[nb - 1] = "-1";
    }

    int[] endPositions = new int[nbSamples];
    String[] endAlleles = new String[nbSamples];

    for (int s = 0; s < nbSamples; s++){
      endPositions[s] = nb; //initialize to last allele
      endAlleles[s] = "*"; //initialize to token value

      for (int i = 0; i < markers.length; i++) {
        Marker m = markers[i];
        String allele = m.getHaplotype(samples[s]);
        String ancestral = m.getAncestral();

        if(noValue(ancestral) || (!allele.equals(ancestral) && !noValue(allele))){//Missing is different
        //if(ancestral.isEmpty() || !allele.equals(ancestral)){//Ignore Missing
          endPositions[s] = i + 1;
          endAlleles[s] = allele;
          break;//stop at the first match
        }
      }
    }
    return new Side(nb, fractions, frequencies, ancestrals, endPositions, endAlleles);
  }

  public static final String S = " ";
  public static final String N = "\n";

  /**
   * Exports an InputFile object to a file
   * @param filename
   * @throws IOException
   */
  public void printToFile(String filename) throws IOException {
    StringBuilder sb = new StringBuilder();
    line(sb, nbSamples, left.getNb(), right.getNb());
    line(sb, left.getFractions());
    line(sb, left.getFrequencies());
    line(sb, right.getFractions());
    line(sb, right.getFrequencies());
    sb.append(mutationRate).append(S).append(mutationModel).append(N);//mix double int
    line(sb, left.getAncestral());
    line(sb, right.getAncestral());
    //Record format:
    //Lp Rp Lv Rv
    // - Lp: 1-based Position of the first Left  Marker with an allele different from the ancestral allele
    // - Rp: 1-based Position of the first Right Marker with an allele different from the ancestral allele
    // - Lv: Allele value at Lp
    // - Rv: Allele value at Rp
    for(int s = 0; s < nbSamples; s++)
      line(sb, left.endPositions[s]+"", right.endPositions[s]+"", left.endAlleles[s], right.endAlleles[s]);

    PrintWriter out = new PrintWriter(new FileWriter(filename));
    out.print(sb);
    out.flush();
    out.close();
  }

  public static void line(StringBuilder sb, int... values){
    sb.append(values[0]);
    for(int i = 1; i < values.length; i++)
      sb.append(S).append(values[i]);
    sb.append(N);
  }

  public static void line(StringBuilder sb, double... values){
    sb.append(values[0]);
    for(int i = 1; i < values.length; i++)
      sb.append(S).append(values[i]);
    sb.append(N);
  }

  public static void line(StringBuilder sb, String... values){
    sb.append(String.join(S, values)).append(N);
  }

  private static boolean isSNP(String s){
    if(s.equals("A") || s.equals("a"))
      return true;
    if(s.equals("C") || s.equals("c"))
      return true;
    if(s.equals("G") || s.equals("g"))
      return true;
    if(s.equals("T") || s.equals("t"))
      return true;
    return false;
  }

  private static boolean isMicrosat(String s){
    try{
      int n = new Integer(s);
      return n > 0;
    } catch(NumberFormatException e) {
      return false;
    }
  }

  /**
   * Class that represents the data on one side of the Marker
   */
  private class Side {
    private final int nb;
    private final double[] fractions;
    private final String[] frequencies;
    private final String[] ancestral;
    private final int[] endPositions;
    private final String[] endAlleles;

    public Side(int nb, double[] fractions, String[] frequencies, String[] ancestral, int[] endPositions, String[] endAlleles) {
      this.nb = nb;
      this.fractions = fractions;
      this.frequencies = frequencies;
      this.ancestral = ancestral;
      this.endPositions = endPositions;
      this.endAlleles = endAlleles;
    }

    public int getNb() {
      return nb;
    }

    public double[] getFractions() {
      return fractions;
    }

    public String[] getFrequencies() {
      return frequencies;
    }

    public String[] getAncestral() {
      return ancestral;
    }

    public int[] getEndPositions() {
      return endPositions;
    }

    public String[] getEndAlleles() {
      return endAlleles;
    }

    public void fromPreinput2Input(){
      for(int i = 0 ; i < nb; i++){
        //if ancestral is -1/end
        if(ancestral[i].equals("-1")) {
          for(int s = 0; s < nbSamples; s++)
            if(endPositions[s] == i+1)
              endAlleles[s] = "1";
        }
        //if ancestral is a snp
        else if(isSNP(ancestral[i])) {
          for(int s = 0; s < nbSamples; s++)
            if(endPositions[s] == i+1) {
              if(ancestral[i].equals(endAlleles[s]))
                endAlleles[s] = "1";
              else if(isSNP(endAlleles[s]))
                endAlleles[s] = "2";
              else
                Message.error("Cannot convert endAllele (position:"+i+", sample:"+s+", allele:"+endAlleles[s]+", ancestral:"+ancestral[i]+")");
            }
          ancestral[i] = "1";
        }
        //if ancestral is a microsat
        else if(isMicrosat(ancestral[i])) {
          int n = new Integer(ancestral[i]);
          int even = n%2;
          for(int s = 0; s < nbSamples; s++){
            if(endPositions[s] == i+1){
              if(ancestral[i].equals(endAlleles[s]))
                endAlleles[s] = "1";
              else if(isMicrosat(endAlleles[s])){
                int v = new Integer(endAlleles[s]);
                if(even != v%2)
                  Message.warning("Parity problem on endAllele (position:"+i+", sample:"+s+", allele:"+endAlleles[s]+", ancestral:"+ancestral[i]+")");
                int d = 1 + Math.abs(n - v) / 2;
                endAlleles[s] = ""+d;
              } else
                Message.error("Cannot convert endAllele (position:"+i+", sample:"+s+", allele:"+endAlleles[s]+", ancestral:"+ancestral[i]+")");
            }
          }
          ancestral[i] = "1";
        }
        //else
        else {
          Message.error("Cannot convert ancestral ["+ancestral[i]+"]");
        }
      }
    }
  }
}
