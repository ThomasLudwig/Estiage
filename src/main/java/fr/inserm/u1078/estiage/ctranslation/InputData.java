package fr.inserm.u1078.estiage.ctranslation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * InputData For Estiage
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-20
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class InputData {
  // Number of markers on each side
  public final int[] nMarkers = new int[2];
  // Number of haplotypes of the longest sequence
  public final int[] nLongest;
  // Maximum number of markers (left of right)
  public final int maxMarkers;
  // Total number of individuals or haplotypes
  public final int nIndividuals;

  /**
   * Recombination fraction for markers on each side from the mutation.
   * Note that the 1st recombination fraction should be set to 0 as it
   * is the recombination fraction between the mutation and itself.
   * Indeed, marker 1 is assumed to be the mutation.
   */
  public final double[][] fractions;// = new double[2][MAX_MARKERS];
  // Mutation rate per generation used in the analysis
  public final double mutationRate;
  // If 1, the stepwise mutation model is used
  public final boolean useStepWiseModel;
  // Position of the 1st marker on each side with a different allele than the ancestral haplotype
  public final int[][] endMarkers;// = new int[2][MAX_INDIVIDUALS];
  // Identity of the allele at the first discordant marker on each side
  public final int[][] endAlleles;// = new int[2][MAX_INDIVIDUALS];
  // Posterior probability of this haplotype reconstruction
  public final double[] postProbability;// = new double[MAX_INDIVIDUALS];
  // Presumed ancestral haplotype on each side, starting from the mutation
  public final int[][] ancestralAlleles;// = new int[2][MAX_MARKERS];
  /**
   * Frequency of the shared allele at each position.
   * Note that freq[k] refers to the frequency of the shared allele at position k-1!
   * Note also that freq[1] and freq[2] are not read and could be set to 0.0.
   */
  public final double[][] frequencies;// = new double[2][MAX_MARKERS];

  // set to false if no data on this side
  public boolean hasLeft;
  // set to false if no data on this side
  public boolean hasRight;

  //* true: alleles are micro sat repetition count, false: alleles are SNP,INDELS
  public boolean isMicrosat = true;

  /**
   * Procedure to read the data from input file
   * July 2001, if markers are available on one side only, the user should
   * put a 0 for the other side
   * January 2003, extension to account for mutations :
   * on each side, the user should specify the alleles at the first unshared
   * markers for each individual
   * April 2006, extension to account for uncertain haplotype : the user
   * should specify the postprob for each haplotype reconstruction
   *
   * @param filename the name of the file to read
   * @throws IOException if the input file can't be read
   * @throws  EstiageException if the input file can't be parsed
   */
  public InputData(String filename) throws IOException, EstiageException {
    boolean tmpLeft = true;
    boolean tmpRight = true;
    String[] line;

    BufferedReader in = new BufferedReader(new FileReader(filename));
    line = in.readLine().split("\\s");

    nIndividuals = Integer.parseInt(line[0]);
    nMarkers[C.LEFT] = Integer.parseInt(line[1]) + 2;
    nMarkers[C.RIGHT] = Integer.parseInt(line[2]) + 2;

    maxMarkers = Math.max(nMarkers[0], nMarkers[1]);
    fractions = new double[2][maxMarkers+1];
    endMarkers = new int[2][nIndividuals];
    endAlleles = new int[2][nIndividuals];
    postProbability = new double[nIndividuals];
    ancestralAlleles = new int[2][maxMarkers+1];
    frequencies = new double[2][maxMarkers+1];

    for (int side : C.SIDES) {
      fractions[side][1] = 0.0;
      line = in.readLine().split("\\s");
      for (int c = 0; c < line.length; c++)
        fractions[side][c + 2] = Double.parseDouble(line[c]);

      line = in.readLine().split("\\s");
      frequencies[side][1] = 0.0;
      frequencies[side][2] = 0.0;
      for (int c = 0; c < line.length; c++)
        frequencies[side][c + 3] = Double.parseDouble(line[c]);
    }
    line = in.readLine().split("\\s");
    mutationRate = Double.parseDouble(line[0]);
    useStepWiseModel = !"0".equals(line[1]);
    line = in.readLine().split("\\s");
    for (int c = 0; c < line.length; c++)
      ancestralAlleles[C.LEFT][c + 1] = readAllele(line[c]);
    line = in.readLine().split("\\s");
    for (int c = 0; c < line.length; c++)
      ancestralAlleles[C.RIGHT][c + 1] = readAllele(line[c]);

    for (int i = 0; i < nIndividuals; i++) {
      line = in.readLine().split("\\s");
      endMarkers[C.LEFT][i] = Integer.parseInt(line[0]) + 1;
      endMarkers[C.RIGHT][i] = Integer.parseInt(line[1]) + 1;
      if (endMarkers[C.LEFT][i] <= 0)
        tmpLeft = false;
      if (endMarkers[C.RIGHT][i] <= 0)
        tmpRight = false;
      endAlleles[C.LEFT][i] = readAllele(line[2]);
      endAlleles[C.RIGHT][i] = readAllele(line[3]);

      if (line.length > 4)
        postProbability[i] = Double.parseDouble(line[4]);
      else
        postProbability[i] = 1.0;
    }
    in.close();
    hasLeft = tmpLeft;
    hasRight = tmpRight;

    this.nLongest = this.findLongest();
  }

  private int readAllele(String s) throws EstiageException {
    if ("-1".equals(s) || s.isEmpty())
      return C.MISSING;
    try {
      int m = Integer.parseInt(s);
      if (!isMicrosat)
        throw new EstiageException("It seems that the input file is mixing microsat and non microsat data");
      return m;
    } catch (NumberFormatException e) {
      this.isMicrosat = false;
      if (useStepWiseModel)
        throw new EstiageException("You cannot use StepWise Model with non microsat data");
      return s.hashCode();
    }
  }

  public int getNMarker(int side) {
    return nMarkers[side];
  }

  public int[] getNLongest() {
    return nLongest;
  }

  /**
   * Procedure to determine the group of haplotypes sharing the longest
   * segment. Left and right sides of the mutation need to be considered
   */
  private int[] findLongest() {
    int[][] count = new int[2][maxMarkers+1];
    int nLongest[] = new int[2];

    for (int side : C.SIDES) {
      for (int j = 1; j <= nMarkers[side]; j++)
        count[side][j] = 0;
      for (int i = 0; i < nIndividuals; i++)
        for (int j = 1; j <= endMarkers[side][i]; j++)
          count[side][j]++;

      int tmp = 0;
      for (int j = 1; j <= nMarkers[side]; j++)
        if (count[side][j] > 1 )
          tmp = j;

      nMarkers[side] = tmp;
      nLongest[side] = count[side][tmp];
      /* Individual markers are recoded */
      for (int i = 0; i < nIndividuals; i++)
        if (endMarkers[side][i] > nMarkers[side])
          endMarkers[side][i] = nMarkers[side];
    }
    return nLongest;
  }

  public int getNIndividuals() {
    return nIndividuals;
  }

  public double[] getFractions(int side) {
    return fractions[side];
  }

  public double getMutationRate() {
    return mutationRate;
  }

  public boolean isUseStepWiseModel() {
    return useStepWiseModel;
  }

  public int getEndMarker(int side, int pos) {
    return endMarkers[side][pos];
  }

  public int getEndAlleles(int side, int pos) {
    return endAlleles[side][pos];
  }

  public double getPostProbability(int pos) {
    return postProbability[pos];
  }

  public int getAncestralAlleles(int side ,int pos) {
    return ancestralAlleles[side][pos];
  }

  public double getFrequencies(int side, int pos) {
    return frequencies[side][pos];
  }
}