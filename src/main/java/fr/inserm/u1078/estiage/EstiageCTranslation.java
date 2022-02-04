package fr.inserm.u1078.estiage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Original EstiAge srouce code
 * Translated from C to Java
 * Some improvement
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-01-14
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class EstiageCTranslation {

  //Constants
  /**
   * Likelihood ratio to stop
   */
  public static final int LIKELIHOOD_STOP_RATIO = 10000;
  /**
   * Maximum number of generations
   */
  public static final int MAX_GENERATIONS = 50000;
  /**
   * index for left side in TAB
   */
  public static final int LEFT = 0;
  /**
   * index for right side in TAB
   */
  public static final int RIGHT = 1;
  /**
   * both sides
   */
  public static final int[] SIDES = {LEFT, RIGHT};

  public static final int MISSING = -1;

  private final InputData input;

  public EstiageCTranslation(String filename) throws IOException, EstiageException {
    input = new InputData(filename);
  }

  public void run(){
    int alFreqOpt = 3;
    int[] nLongest = findLongest(false);
    /* Final Maximum Log Likelihood */
    try {
      Results res = findMaxLike(alFreqOpt, nLongest);
      res.print();
    } catch (EstiageException e) {
      System.err.println(e.getMessage());
    }
  }

  private double S(int x, int nGenerations, double[] vector) {
    if (x > 0) {
      double recFrac = vector[x];
      return Math.pow(1.0 - recFrac, nGenerations);
    }
    return 0;
  }

  private double F(int x, int nGenerations, double[] vector) {
    if (x > 1)
      return S(x - 1, nGenerations, vector) - S(x, nGenerations, vector);
    return 0;
  }

  /**
   * Compute n!
   * @param n an integer
   * @return factorial(n)
   */
  private double fact(int n) {
    double prod = 1.0;
    //fact(0) = fact(1) = 1;
    for (int i = 2; i <= n; i++)
      prod *= i;
    return prod;
  }

  /**
   * Probability that there was no mutation after several generation
   * @param x
   * @param nGeneration the number of generation
   * @return
   */
  private double U(int x, int nGeneration) {
    return Math.pow(1.0 - input.mutationRate, nGeneration * (x - 1));
  }

  private double H(int x, int nGenerations, double[] vector, int a1, int a2) {
    double pmut = 0;

    if ((a2 >= 0) && (a1 >= 0)) {
      int obsstep = Math.abs(a2 - a1);
      if (input.useStepWiseModel) {
        double lambda = input.mutationRate * (double) nGenerations;
        if (lambda > 0)
          pmut = Math.pow(lambda, obsstep) * Math.exp(-lambda) / fact(obsstep);
        else
          pmut = 0.0;
      } else
        pmut = 1.0 - Math.pow(1.0 - input.mutationRate, nGenerations);
    }
    return F(x, nGenerations, vector) + pmut * S(x, nGenerations, vector);
  }

  /**
   * Computation of the likelihood over the whole sample.
   *
   * @param side        LEFT:0 RIGHT:1
   * @param nGeneration number of Generation
   * @param ny          number of individuals belonging to the G1 group
   */
  private double[] like(int side, int nGeneration, int ny) {
    double lik1b = 0, lik2b = 0, lik3b = 0;
    /* Computations for individuals belonging to the G1 group */
    /* First case: one individual carries the ancestral haplotype */
    for (int i = 0; i < input.nIndividuals; i++) {
      /* We need to determine who is this individual */
      double likTmp1 = 1.0;
      int xi = input.endMarkers[side][i];
      if (xi == input.nMarkers[side]) {
        double likTmp2 = 0;
        double likTmp3 = 0;
        double contribAnc = U(xi, nGeneration) * S(xi, nGeneration, input.rec[side]);
        contribAnc = Math.pow(contribAnc, input.postProbability[i]);
        double totLikTmp2 = 0.0;
        double totLikTmp3 = 0.0;
        double totTotLikTmp4 = 0.0;
        for (int j = 0; j < input.nIndividuals; j++) {
          double totLikTmp4 = 0.0;
          int xj = input.endMarkers[side][j];
          if ((xj == input.nMarkers[side]) && (j != i)) {
            double p1 = (xj > 2) ? input.freq[side][xj]: 0.0;
            double p2 = (xj > 2) ? input.freq[side][xj - 1] : 0.0;
            double contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], input.endAlleles[side][i], input.endAlleles[side][j]);
            contrib = Math.pow(contrib, input.postProbability[j]);
            likTmp1 *= contrib;
            if (ny >= 2) {
              likTmp2 = U(xj - 1, nGeneration) * p1 * F(xj - 1, nGeneration, input.rec[side]);
              likTmp2 = Math.pow(likTmp2, input.postProbability[j]);
              likTmp3 = U(xj - 1, nGeneration) * p1 * p2 * F(xj - 2, nGeneration, input.rec[side]);
              likTmp3 = Math.pow(likTmp3, input.postProbability[j]);
              for (int k = 0; k < input.nIndividuals; k++)
                if (input.endMarkers[side][k] == input.nMarkers[side] && (k != i) && (k != j)) {
                  contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], input.endAlleles[side][i], input.endAlleles[side][k]);
                  contrib = Math.pow(contrib, input.postProbability[k]);
                  likTmp2 *= contrib;
                  likTmp3 *= contrib;
                }
            }
            if (ny >= 3)
              for (int l = j + 1; l < input.nIndividuals; l++)
                if ((input.endMarkers[side][l] == input.nMarkers[side]) && (l != i)) {
                  double likTmp4 = Math.pow(U(xj - 1, nGeneration) * p1 * F(xj - 1, nGeneration, input.rec[side]), 2);
                  likTmp4 = Math.pow(likTmp4, input.postProbability[l] * input.postProbability[j]);
                  for (int m = 0; m < input.nIndividuals; m++)
                    if (input.endMarkers[side][m] == input.nMarkers[side] && m != l && (m != i) && (m != j)) {
                      contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], input.endAlleles[side][i], input.endAlleles[side][m]);
                      contrib = Math.pow(contrib, input.postProbability[m]);
                      likTmp4 *= contrib;
                    }
                  totLikTmp4 += likTmp4;
                }
            totLikTmp2 += likTmp2;
            totLikTmp3 += likTmp3;
            totTotLikTmp4 += totLikTmp4;
          }
        }  /* Loop on individual j */
        lik1b += contribAnc * likTmp1;
        lik2b += contribAnc * (likTmp1 + totLikTmp2);
        lik3b += contribAnc * (likTmp1 + totLikTmp2 + totLikTmp3 + totTotLikTmp4);
      }
    } /* Loop on individual i, the ancestor */
    /* Computations for individuals belonging to the G1 group */
    /* Second case: no individual carries the ancestral haplotype */
    double likTmp1 = 1.0;
    double totLikTmp2 = 0.0;
    double totLikTmp3 = 0.0;
    double totTotLikTmp4 = 0.0;
    for (int j = 0; j < input.nIndividuals; j++) {
      double totLikTmp4 = 0.0;
      int xj = input.endMarkers[side][j];
      if (xj == input.nMarkers[side]) {
        double likTmp2 = 0;
        double likTmp3 = 0;
        double p1 = (xj > 2) ? input.freq[side][xj]: 0.0;
        double p2 = (xj > 2) ? input.freq[side][xj - 1] : 0.0;
        double contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], MISSING, input.endAlleles[side][j]);
        contrib = Math.pow(contrib, input.postProbability[j]);
        likTmp1 *= contrib;
        if (ny >= 2) {
          likTmp2 = U(xj - 1, nGeneration) * p1 * F(xj - 1, nGeneration, input.rec[side]);
          likTmp2 = Math.pow(likTmp2, input.postProbability[j]);
          likTmp3 = U(xj - 1, nGeneration) * p1 * p2 * F(xj - 2, nGeneration, input.rec[side]);
          likTmp3 = Math.pow(likTmp3, input.postProbability[j]);
          for (int k = 0; k < input.nIndividuals; k++)
            if (input.endMarkers[side][k] == input.nMarkers[side] && k != j) {
              contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], MISSING, input.endAlleles[side][k]);
              contrib = Math.pow(contrib, input.postProbability[k]);
              likTmp2 *= contrib;
              likTmp3 *= contrib;
            }
        }
        if (ny >= 3)
          for (int l = j + 1; l < input.nIndividuals; l++)
            if ((input.endMarkers[side][l] == input.nMarkers[side]) && (l != input.nIndividuals-1)) {//l != i-1, here i is always nind, from the end of top level loop
              double likTmp4 = Math.pow(U(xj - 1, nGeneration) * p1 * F(xj - 1, nGeneration, input.rec[side]), 2);
              likTmp4 = Math.pow(likTmp4, input.postProbability[l]);
              for (int m = 0; m < input.nIndividuals; m++)
                if (input.endMarkers[side][m] == input.nMarkers[side] && m != l && (m != input.nIndividuals-1) && (m != j)) { //m != i-1, here i is always nind, from the end of top level loop
                  contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], MISSING, input.endAlleles[side][m]);
                  contrib = Math.pow(contrib, input.postProbability[m]);
                  likTmp4 *= contrib;
                }
              totLikTmp4 += likTmp4;
            }

        totLikTmp2 += likTmp2;
        totLikTmp3 += likTmp3;
        totTotLikTmp4 += totLikTmp4;
      }
    }  /* Loop on individual j */
    double lik1a = likTmp1;
    double lik2a = likTmp1 + totLikTmp2;
    double lik3a = likTmp1 + totLikTmp2 + totLikTmp3 + totTotLikTmp4;

    /* Computations for individuals belonging to the G2 group */
    double lik1c = 1.0;
    double lik2c = 1.0;
    double lik3c = 1.0;
    for (int j = 0; j < input.nIndividuals; j++) {
      int xj = input.endMarkers[side][j];
      if (xj != input.nMarkers[side]) {
        double p1 = (xj > 2) ? input.freq[side][xj]: 0.0;
        double p2 = (xj >= 3) ? input.freq[side][xj - 1] : 0.0;
        double contrib = U(xj - 1, nGeneration) * H(xj, nGeneration, input.rec[side], input.ancestralAlleles[side][xj], input.endAlleles[side][j]);
        lik1c *= Math.pow(contrib, input.postProbability[j]);
        lik2c *= Math.pow(contrib + p1 * F(xj - 1, nGeneration, input.rec[side]), input.postProbability[j]);
        lik3c *= Math.pow(contrib + p1 * F(xj - 1, nGeneration, input.rec[side]) + p1 * p2 * F(xj - 2, nGeneration, input.rec[side]), input.postProbability[j]);
      }
    }
    return new double[]{
            0,
            (lik1a + lik1b) * lik1c,
            (lik2a + lik2b) * lik2c,
            (lik3a + lik3b) * lik3c
    };
  }

  /**
   * Procedure to determine the group of haplotypes sharing the longest
   * segment. Left and right sides of the mutation need to be considered
   *
   * @param byGroup do it by group ?
   * @return Number of subjects with the longest haplotype on each side
   */
  private int[] findLongest(boolean byGroup) {
    int[][] count = new int[2][input.maxMarkers+1];
    int nLongest[] = new int[2];
    nLongest[LEFT] = 0;
    nLongest[RIGHT] = 0;

    for (int side : SIDES)
      for (int j = 1; j <= input.nMarkers[side]; j++)
        count[side][j] = 0;

    for (int i = 0; i < input.nIndividuals; i++)
      if (!byGroup)
        for (int side : SIDES)
          for (int j = 1; j <= input.endMarkers[side][i]; j++)
            count[side][j]++;

    /* Determination of nmar[side] based on the data */
    for (int side : SIDES) {
      int tmp = 0;

      for (int j = 1; j <= input.nMarkers[side]; j++)
        if ((count[side][j] > 1) && (j > tmp))
          tmp = j;

      input.nMarkers[side] = tmp;
      nLongest[side] = count[side][tmp];
      /* Individual markers are recoded */
      for (int i = 0; i < input.nIndividuals; i++)
        if (!byGroup)
          if (input.endMarkers[side][i] > input.nMarkers[side])
            input.endMarkers[side][i] = input.nMarkers[side];
    }
    return nLongest;
  }

  /**
   * Computation for subjects sharing the longest haplotype
   *
   * @param ng
   * @param typelik
   * @return
   */
  private double totLike(int ng, int typelik, final int[] nLongest) {
    double likeTot = 1;
    if (input.hasLeft) {
      double[] likeCalc = like(LEFT, ng, nLongest[LEFT]);
      likeTot = likeCalc[typelik];
    }

    if (input.hasRight) {
      double[] likeCalc = like(RIGHT, ng, nLongest[RIGHT]);
      likeTot *= likeCalc[typelik];
    }
    return likeTot;
  }

  /**
   * Main computation
   *
   * @param typelik alFreqOpt
   * @return results of Estiage
   * @throws EstiageException
   */
  private Results findMaxLike(int typelik, final int[] nLongest) throws EstiageException {
    /** Number of generations from common ancestor */
    int nGen = 1;
    /** Sum of the likelihood over nGen */
    double pTot = 0;
    /** Temporary variable to compute 95% CI */
    double pTic = 0;
    /**
     * Keep the maximum likelihood at each incrementation of nGen
     */
    double pMax = totLike(nGen, typelik, nLongest);
    /**
     * Array to keep the likelihood of n
     */
    double[] pGenTot = new double[MAX_GENERATIONS + 1];
    int nMax = nGen;
    for (int gen = 1; gen <= MAX_GENERATIONS; gen++)
      pGenTot[gen] = 0.0;
    /* Computation of the maximum likelihood */
    pGenTot[nGen] = totLike(nGen, typelik, nLongest);
    while (((pMax / pGenTot[nGen]) < LIKELIHOOD_STOP_RATIO) && (nGen <= MAX_GENERATIONS)) {
      if (pGenTot[nGen] >= pMax) {
        pMax = pGenTot[nGen];
        nMax = nGen;
      }
      nGen++;
      pGenTot[nGen] = totLike(nGen, typelik, nLongest);
    }

    if (nGen > MAX_GENERATIONS) {
      String message = "Maximum number of iterations reached";
      for (int i = 0; i < input.nIndividuals; i++)
        message += "\n" + i + " " + input.endMarkers[LEFT][i] + " " + input.endMarkers[RIGHT][i];
      throw new EstiageException(message);
    }

    int nEnd = nGen;
    int nInf = 0;
    int nSup = 0;
    for (int gen = 1; gen <= nEnd; gen++)
      pTot += pGenTot[gen];
    for (int gen = 1; gen <= nEnd; gen++) {
      pTic += pGenTot[gen] / pTot;
      if ((pTic > 0.025) && (nInf <= 0))
        nInf = gen;
      if (pTic > 0.975) {
        nSup = gen;
        break;
      }
    }
    double lnTotLike = Math.log(pMax);
    ;
    return new Results(nMax, nEnd, nInf, nSup, lnTotLike);
  }


  private static class InputData {
    /**
     * Number of markers on each side
     */
    private final int[] nMarkers = new int[2];
    /**
     * Maximum number of markers (left of right)
     */
    private final int maxMarkers;
    /**
     * Total number of individuals or haplotypes
     */
    private final int nIndividuals;

    /**
     * Recombination fraction for markers on each side from the mutation.
     * Note that the 1st recombination fraction should be set to 0 as it
     * is the recombination fraction between the mutation and the mutation.
     * Indeed, marker 1 is assumed to be the mutation.
     */
    private final double[][] rec;// = new double[2][MAX_MARKERS];
    /**
     * Mutation rate per generation used in the analysis
     */
    private final double mutationRate;
    /**
     * If 1, the stepwise mutation model is used
     */
    private final boolean useStepWiseModel;
    /**
     * Position of the 1st marker on each side with a different allele than the ancestral haplotype
     */
    private final int[][] endMarkers;// = new int[2][MAX_INDIVIDUALS];
    /**
     * Identity of the allele at the first discordant marker on each side
     */
    private final int[][] endAlleles;// = new int[2][MAX_INDIVIDUALS];
    /**
     * Posterior probability of this haplotype reconstruction
     */
    private final double[] postProbability;// = new double[MAX_INDIVIDUALS];
    /**
     * Presumed ancestral haplotype on each side, starting from the mutation
     */
    private final int[][] ancestralAlleles;// = new int[2][MAX_MARKERS];
    /**
     * Frequency of the shared allele at each position.
     * Note that freq[k] refers to the frequency of the shared allele at position k-1!
     * Note also that freq[1] and freq[2] are not read and could be set to 0.0.
     */
    private final double[][] freq;// = new double[2][MAX_MARKERS];

    //Various
    /**
     * set to false if no data on this side
     */
    boolean hasLeft = true;
    /**
     * set to false if no data on this side
     */
    boolean hasRight = true;

    /**
     * true: alleles are micro sat repetition count, false: alleles are SNP,INDELS
     */
    boolean isMicrosat = true;

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
     * @param filename the file to read
     * @throws IOException
     */
    InputData(String filename) throws IOException, EstiageException {
      String[] line;

      BufferedReader in = new BufferedReader(new FileReader(filename));
      line = in.readLine().split("\\s");

      nIndividuals = Integer.parseInt(line[0]);
      nMarkers[LEFT] = Integer.parseInt(line[1]) + 2;
      nMarkers[RIGHT] = Integer.parseInt(line[2]) + 2;

      maxMarkers = Math.max(nMarkers[0], nMarkers[1]);
      rec = new double[2][maxMarkers+1];
      endMarkers = new int[2][nIndividuals];
      endAlleles = new int[2][nIndividuals];
      postProbability = new double[nIndividuals];
      ancestralAlleles = new int[2][maxMarkers+1];
      freq = new double[2][maxMarkers+1];

      for (int side : SIDES) {
        rec[side][1] = 0.0;
        line = in.readLine().split("\\s");
        for (int c = 0; c < line.length; c++)
          rec[side][c + 2] = Double.parseDouble(line[c]);

        line = in.readLine().split("\\s");
        freq[side][1] = 0.0;
        freq[side][2] = 0.0;
        for (int c = 0; c < line.length; c++)
          freq[side][c + 3] = Double.parseDouble(line[c]);
      }
      line = in.readLine().split("\\s");
      mutationRate = Double.parseDouble(line[0]);
      useStepWiseModel = !"0".equals(line[1]);
      line = in.readLine().split("\\s");
      for (int c = 0; c < line.length; c++)
        ancestralAlleles[LEFT][c + 1] = readAllele(line[c]);
      line = in.readLine().split("\\s");
      for (int c = 0; c < line.length; c++)
        ancestralAlleles[RIGHT][c + 1] = readAllele(line[c]);

      for (int i = 0; i < nIndividuals; i++) {
        line = in.readLine().split("\\s");
        endMarkers[LEFT][i] = Integer.parseInt(line[0]) + 1;
        endMarkers[RIGHT][i] = Integer.parseInt(line[1]) + 1;
        if (endMarkers[LEFT][i] <= 0)
          hasLeft = false;
        if (endMarkers[RIGHT][i] <= 0)
          hasRight = false;
        endAlleles[LEFT][i] = readAllele(line[2]);
        endAlleles[RIGHT][i] = readAllele(line[3]);

        if (line.length > 4)
          postProbability[i] = Double.parseDouble(line[4]);
        else
          postProbability[i] = 1.0;
      }
      in.close();
    }

    private int readAllele(String s) throws EstiageException {
      if ("-1".equals(s) || s.isEmpty())
        return MISSING;
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
  }

  /**
   * Class representing the Results of Estiage.
   */
  public static class Results {
    /* Number of generations:maxlike,end of process,95% CI smaller limit, 95% CI upper limit */
    private final int max, end, inf, sup;
    private final double lnTotLike;

    public Results(int max, int end, int inf, int sup, double lntotlike) {
      this.max = max;
      this.end = end;
      this.inf = inf;
      this.sup = sup;
      this.lnTotLike = lntotlike;
    }

    /**
     * Prints results on StdOut
     */
    public void print() {
      System.out.println("n = " + max + ", " +
              "nend = " + end + ", " +
              "ninf = " + inf + ", " +
              "nsup = " + sup + ", " +
              "likelihood = " + lnTotLike);
    }
  }

  /**
   * Exception during the call to Estiage
   */
  public static class EstiageException extends Exception {
    public EstiageException(String message) {
      super(message);
    }
  }
}

