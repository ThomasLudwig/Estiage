package fr.inserm.u1078.estiage.ctranslation;

import fr.inserm.u1078.estiage.MathLib;

/**
 * Constants
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-20
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class C {
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

  /**
   * MISSING ALLELE
   */
  public static final int MISSING = -1;

  /**
   *
   * @param rank rank of the last marker
   * @param nGenerations the number of generation
   * @param fractions the recombination fractions
   * @return 0 if the rank is negative or null, (1-recombinationFraction)^nGenerations otherwise
   */
  public static double S(int rank, int nGenerations, final double[] fractions) {
    return rank <= 0 ? 0 : Math.pow(1.0 - fractions[rank], nGenerations);
  }

  /**
   *
   * @param rank rank of the last marker
   * @param nGenerations the number of generation
   * @param fractions the recombination fractions
   * @return 0 if rank <= 1, S(rank-1) - S(rank) otherwise
   */
  public static double F(int rank, int nGenerations, final double[] fractions) {
    return rank <= 1 ? 0 : S(rank - 1, nGenerations, fractions) - S(rank, nGenerations, fractions);
  }

  /**
   * Probability that there was no mutation after several generation
   * @param mutationRate the mutation Rate in the model
   * @param rank rank of the last marker
   * @param nGenerations the number of generation
   * @return
   */
  public static double U(double mutationRate, int rank, int nGenerations) {
    return Math.pow(1.0 - mutationRate, nGenerations * (rank - 1));
  }

  /**
   *
   * @param mutationRate the mutation Rate in the model
   * @param rank rank of the last marker
   * @param nGenerations the number of generation
   * @param fractions the recombination fractions
   * @param endAllele1 the first end allele
   * @param endAllele2 the second end allele
   * @param stepWise true, if the stepwise model is used
   * @return
   */
  public static double H(double mutationRate, int rank, int nGenerations, final double[] fractions, int endAllele1, int endAllele2, boolean stepWise) {
    double pMut;

    if ((endAllele2 >= 0) && (endAllele1 >= 0)) {
      if (stepWise) {
        double lambda = mutationRate * nGenerations;
        if (lambda > 0) {
          int obsStep = Math.abs(endAllele2 - endAllele1);
          pMut = Math.pow(lambda, obsStep) * Math.exp(-lambda) / MathLib.fact(obsStep);
        }
        else
          pMut = 0.0;
      } else
        pMut = 1.0 - Math.pow(1.0 - mutationRate, nGenerations);
    } else
      pMut = 0.0;
    return F(rank, nGenerations, fractions) + pMut * S(rank, nGenerations, fractions);
  }
}
