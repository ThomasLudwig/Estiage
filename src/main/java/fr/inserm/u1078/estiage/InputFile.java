package fr.inserm.u1078.estiage;

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

  public static final int LEFT = 0;
  public static final int RIGHT = 1;
  public static final int IP = 0;
  public static final int IV = 2;

  public static boolean noValue(String s){
    if("-1".equals(s))
      return true;
    return s.isEmpty();
  }

  private Side load(/*int side, */Marker[] markers, String[] samples){//}, double[] fractions, double[] frequencies, String[] ancestrals){
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
  }
}
