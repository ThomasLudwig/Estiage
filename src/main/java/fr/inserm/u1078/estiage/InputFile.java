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
  public static final String SEP = " ";
  public static final String NL = "\n";
  private final int nbSamples;
  private final double mutationRate;
  private final int mutationModel;

  private final Side left;
  private final Side right;

  /**
   * Creates an InputFile object
   * @param raw the name of the source raw file
   * @param mutationModel the mutation model to use
   * @param mutationRate the mutation rate
   */
  public InputFile(TSVFile raw, int mutationModel, double mutationRate) {
    this.nbSamples = raw.getSamples().length;
    this.mutationRate = mutationRate;
    this.mutationModel = mutationModel;

    this.left = Side.build(raw.getLeftMarkers(), raw.getSamples());
    this.right = Side.build(raw.getRightMarkers(), raw.getSamples());
  }

  /**
   * Constructor. Builds an InputFile object from a serialized input file
   * @param preInputFile the preInput File
   * @throws EstiageFormatException if the PreInput File can't be parsed
   * @throws IOException if the preInputFile can't be read
   */
  public InputFile(String preInputFile) throws  EstiageFormatException, IOException {
    UniversalReader in = new UniversalReader(preInputFile);
    String[] f;
    int line = 0;
    try {
      //header
      line++;
      f = in.readLine().split(SEP, -1);
      this.nbSamples = Integer.parseInt(f[0]);
      int nbLeft = Integer.parseInt(f[1]);
      int nbRight = Integer.parseInt(f[2]);

      int[] positionsLeft = new int[nbSamples];
      String[] endLeft = new String[nbSamples];
      int[] positionsRight = new int[nbSamples];
      String[] endRight = new String[nbSamples];

      //fractions left
      line++;
      double[] fractionsLeft = parseDouble(in.readLine(), nbLeft);

      //frequencies left
      line++;
      String[] frequenciesLeft = splitString(in.readLine(), nbLeft);

      //fractions right
      line++;
      double[] fractionsRight = parseDouble(in.readLine(), nbRight);

      //frequencies right
      line++;
      String[] frequenciesRight = splitString(in.readLine(), nbRight);

      //mutation
      line++;
      f = in.readLine().split(SEP, -1);
      this.mutationRate = Double.parseDouble(f[0]);
      this.mutationModel = Integer.parseInt(f[1]);

      //ancestral left
      line++;
      String[] ancestralLeft = splitString(in.readLine(), nbLeft);

      //ancestral right
      line++;
      String[] ancestralRight = splitString(in.readLine(), nbRight);

      // samples
      for(int i = 0 ; i < nbSamples; i++) {
        line++;
        f = in.readLine().split(SEP, -1);
        positionsLeft[i] = Integer.parseInt(f[0]);
        positionsRight[i] = Integer.parseInt(f[1]);
        endLeft[i] = f[2];
        endRight[i] = f[3];
      }

      this.left = new Side(nbLeft, fractionsLeft, frequenciesLeft, ancestralLeft, positionsLeft, endLeft);
      this.right = new Side(nbRight, fractionsRight, frequenciesRight, ancestralRight, positionsRight, endRight);
    } catch(NullPointerException | NumberFormatException | ArrayIndexOutOfBoundsException e) {
      throw new EstiageFormatException("Unable to parse input file ["+preInputFile+"] on line ["+line+"]", e);
    }
    in.close();
  }

  /**
   * Splits a line into an array of double
   * @param line the input line
   * @param n the number of expected values
   * @return the values
   */
  private static double[] parseDouble(String line, int n) {
    double[] ret = new double[n];
    String[] f = line.split(SEP, -1);
    for(int i = 0 ; i < n; i++)
      ret[i] = Double.parseDouble(f[i]);
    return ret;
  }

  /**
   * Splits a line into an array of string
   * @param line the input line
   * @param n the number of expected values
   * @return the values
   */
  private static String[] splitString(String line, int n) {
    String[] ret = new String[n];
    String[] f = line.split(SEP, -1);
    System.arraycopy(f, 0, ret, 0, n);
    return ret;
  }

  /**
   * Converts a preInput file to an input file
   */
  public void fromPreInput2Input(){
    left.fromPreInput2Input();
    right.fromPreInput2Input();
  }

  /**
   * Check if a string contains no values
   * @param s the string to check
   * @return true, if the string == "-1" or ""
   */
  public static boolean noValue(String s){
    if("-1".equals(s))
      return true;
    return s.isEmpty();
  }

  /**
   * Exports an InputFile object to a file
   * @param filename the name of the file
   * @throws IOException if the file can't be written
   */
  public void export(String filename) throws IOException {
    StringBuilder sb = new StringBuilder();
    append(sb, nbSamples, left.getNb(), right.getNb());
    append(sb, left.getFractions());
    append(sb, left.getFrequencies());
    append(sb, right.getFractions());
    append(sb, right.getFrequencies());
    sb.append(mutationRate).append(SEP).append(mutationModel).append(NL);//mix double int
    append(sb, left.getAncestral());
    append(sb, right.getAncestral());
    //Record format:
    //Lp Rp Lv Rv
    // - Lp: 1-based Position of the first Left  Marker with an allele different from the ancestral allele
    // - Rp: 1-based Position of the first Right Marker with an allele different from the ancestral allele
    // - Lv: Allele value at Lp
    // - Rv: Allele value at Rp
    for(int s = 0; s < nbSamples; s++)
      append(sb, left.endPositions[s]+"", right.endPositions[s]+"", left.endAlleles[s], right.endAlleles[s]);

    PrintWriter out = new PrintWriter(new FileWriter(filename));
    out.print(sb);
    out.flush();
    out.close();
  }

  /**
   * append a series of integer values to a line (in StringBuilder format)
   * @param sb the line in StringBuilderFormat
   * @param values the series of integer values
   */
  public static void append(StringBuilder sb, int... values){
    sb.append(values[0]);
    for(int i = 1; i < values.length; i++)
      sb.append(SEP).append(values[i]);
    sb.append(NL);
  }

  /**
   * append a series of double values to a line (in StringBuilder format)
   * @param sb the line in StringBuilderFormat
   * @param values the series of double values
   */
  public static void append(StringBuilder sb, double... values){
    sb.append(values[0]);
    for(int i = 1; i < values.length; i++)
      sb.append(SEP).append(values[i]);
    sb.append(NL);
  }

  /**
   * append a series of String values to a line (in StringBuilder format)
   * @param sb the line in StringBuilderFormat
   * @param values the series of String values
   */
  public static void append(StringBuilder sb, String... values){
    sb.append(String.join(SEP, values)).append(NL);
  }

  /**
   * checks if a string value is a Sequence Variant (ie not a number)
   * @param s the value
   * @return if the value is compatible with a SNP
   */
  private static boolean isSequenceVariant(String s){
    final int token = -589798156;
    int value = token;
    try{
      value = Integer.parseInt(s);
    } catch (NumberFormatException ignore){
      //ignore
    }
    return value == token;
  }

  /**
   * Checks if a string is a microsat value (ie a positive integer)
   * @param s the value
   * @return true if value is compatible with microsat
   */
  private static boolean isMicrosatVariant(String s){
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
  private static class Side {
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

    /**
     * Creates a Side object from Markers and Samples
     * @param markers an array of markers
     * @param samples an array of samples
     * @return the Side object
     */
    public static Side build(Marker[] markers, String[] samples){
      int M = markers.length;
      final int S = samples.length;
      boolean ambiguous = !noValue(markers[M-1].getAncestral());

      if(ambiguous)
        M++;

      double[] fractions = new double[M];
      String[] frequencies = new String[M];
      String[] ancestrals = new String[M];

      for (int i = 0; i < markers.length; i++) {
        Marker m = markers[i];
        String ancestral = m.getAncestral();
        fractions[i] = m.getRecombinationFraction();
        frequencies[i] = ""+m.getFrequency();
        ancestrals[i] = ancestral;
      }

      if(ambiguous){//If the last marker has a clear ancestral allele, add another marker with a fraction of .5 after this one
        fractions[M - 1] = .5;
        frequencies[M - 1] = "";
        ancestrals[M - 1] = "-1";
      }

      int[] endPositions = new int[S];
      String[] endAlleles = new String[S];

      for (int s = 0; s < S; s++){
        endPositions[s] = M; //initialize to last allele
        endAlleles[s] = "*"; //initialize to token value

        for (int i = 0; i < markers.length; i++) {
          Marker m = markers[i];
          String allele = m.getAllele(samples[s]);
          String ancestral = m.getAncestral();

          if(noValue(ancestral) || (!allele.equals(ancestral) && !noValue(allele))){//Missing is different
            //if(ancestral.isEmpty() || !allele.equals(ancestral)){//Ignore Missing
            endPositions[s] = i + 1;
            endAlleles[s] = allele;
            break;//stop at the first match
          }
        }
      }
      return new Side(M, fractions, frequencies, ancestrals, endPositions, endAlleles);
    }

    /**
     * Gets the number of markers
     * @return the number of markers
     */
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

    /**
     * Transform the PreInput Side to an InputSide
     */
    public void fromPreInput2Input(){
      final int nbSamples = endAlleles.length;
      for(int i = 0 ; i < nb; i++){
        //if ancestral is -1/end
        if(ancestral[i].equals("-1")) {
          for(int s = 0; s < nbSamples; s++)
            if(endPositions[s] == i+1)
              endAlleles[s] = "1";
        }
        //if ancestral is a sequence variant
        else if(isSequenceVariant(ancestral[i])) {
          for(int s = 0; s < nbSamples; s++)
            if(endPositions[s] == i+1) {
              if(ancestral[i].equals(endAlleles[s]))
                endAlleles[s] = "1";
              else if(isSequenceVariant(endAlleles[s]))
                endAlleles[s] = "2";
              else
                Message.error("Cannot convert endAllele (position:"+i+", sample:"+s+", allele:"+endAlleles[s]+", ancestral:"+ancestral[i]+")");
            }
          ancestral[i] = "1";
        }
        //if ancestral is a microsat
        else if(isMicrosatVariant(ancestral[i])) {
          int n = new Integer(ancestral[i]);
          int even = n%2;
          for(int s = 0; s < nbSamples; s++){
            if(endPositions[s] == i+1){
              if(ancestral[i].equals(endAlleles[s]))
                endAlleles[s] = "1";
              else if(isMicrosatVariant(endAlleles[s])){
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