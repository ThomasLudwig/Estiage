package fr.inserm.u1078.estiage.ctranslation;

import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.IOException;

/**
 * Original EstiAge source code
 * Translated from C to Java
 * Some improvement
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-01-14
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Estiage {

  /**
   * Runs the EstiAge algorithm on the input file
   * @param filename the name of the input file
   * @throws IOException if the file can't be read
   * @throws EstiageException if the file can't be parsed
   */
  public static void run(String filename) throws IOException, EstiageException {
    final InputData input = new InputData(filename);
    int[] nLongest = input.getNLongest();
    /* Final Maximum Log Likelihood */
    try {
      Results res = findMaxLike(input, nLongest);
      res.print();
    } catch (EstiageException e) {
      System.err.println(e.getMessage());
    }
  }

  /**
   * Computation of the likelihood over the whole sample.
   * @param input       the input data
   * @param side        LEFT:0 RIGHT:1
   * @param nGeneration number of Generation
   * @param n           number of individuals belonging to the G1 group
   * @return the likelihood
   */
  private static double like(InputData input, int side, int nGeneration, int n) {
    /* no individual carries the ancestral haplotype */
    final double likA = likeA(input, side, nGeneration, n);
    /* Computations for individuals belonging to the G1 group (carriers) */
    final double likB = likeB(input, side, nGeneration, n);
    /* Computations for individuals belonging to the G2 group (non carriers) */
    final double likC = likeC(input, side, nGeneration);
    Message.debug("a="+likA);
    Message.debug("b="+likB);
    Message.debug("c="+likC);
    Message.debug("like="+(likA + likB) * likC);
    return (likA + likB) * likC;
  }

  /**
   * No individual Carries the ancestral haplotype
   * @param input       the input data
   * @param side        LEFT:0 RIGHT:1
   * @param nGeneration number of Generation
   * @param n          number of individuals belonging to the G1 group
   * @return the likelihood
   */
  private static double likeA(InputData input, int side, int nGeneration, int n) {
    /* Second case: no individual carries the ancestral haplotype */
    int nMarker = input.getNMarker(side);
    double likTmp1 = 1.0;
    double totLikTmp2 = 0.0;
    double totLikTmp3 = 0.0;
    double totTotLikTmp4 = 0.0;
    /* Loop on individual j */
    for (int j = 0; j < input.getNIndividuals(); j++) {
      double totLikTmp4 = 0.0;
      final int endMarkerJ = input.getEndMarker(side, j);
      if (endMarkerJ == nMarker) {
        double likTmp2 = 0;
        double likTmp3 = 0;
        final double p1 = (endMarkerJ > 2) ? input.getFrequencies(side, endMarkerJ) : 0.0;
        final double p2 = (endMarkerJ > 2) ? input.getFrequencies(side, endMarkerJ - 1) : 0.0;
        final double u = C.U(input.getMutationRate(), endMarkerJ - 1, nGeneration);
        final double hj = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), C.MISSING, input.getEndAlleles(side, j), input.isUseStepWiseModel());
        likTmp1 *= Math.pow(u * hj, input.getPostProbability(j));

        //If more than 1 G1
        if (n > 1) {
          final double f1 = C.F(endMarkerJ - 1, nGeneration, input.getFractions(side));
          final double f2 = C.F(endMarkerJ - 2, nGeneration, input.getFractions(side));
          //only last value wil be kept
          likTmp2 = Math.pow(u * p1 * f1, input.getPostProbability(j));
          likTmp3 = Math.pow(u * p1 * p2 * f2, input.getPostProbability(j));
          Message.debug("j=" + j + " t2=" + likTmp2 + " t3=" + likTmp3);
          for (int k = 0; k < input.getNIndividuals(); k++) {
            if (input.getEndMarker(side, k) == nMarker && k != j) {
              final double hk = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), C.MISSING, input.getEndAlleles(side, k), input.isUseStepWiseModel());
              final double contrib = Math.pow(u * hk, input.getPostProbability(k));
              likTmp2 *= contrib;
              likTmp3 *= contrib;
              Message.debug("j=" + j + " k=" + k + " t2=" + likTmp2 + " t3=" + likTmp3);
            }
          }


          //If more than 2 G1
          if (n > 2) {
            for (int l = j + 1; l < input.getNIndividuals(); l++) {
              if ((input.getEndAlleles(side, l) == nMarker) && (l != input.getNIndividuals() - 1)) {//l != i-1, here i is always nind, from the end of top level loop
                final double tmpCF = Math.pow(u * p1 * f1, 2);
                double likTmp4 = Math.pow(tmpCF, input.getPostProbability(l));
                for (int m = 0; m < input.getNIndividuals(); m++)
                  if (input.getEndAlleles(side, m) == nMarker && (m != l) && (m != input.getNIndividuals() - 1) && (m != j)) { //m != i-1, here i is always nind, from the end of top level loop
                    final double hm = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), C.MISSING, input.getEndAlleles(side, m), input.isUseStepWiseModel());
                    final double contrib = Math.pow(u * hm, input.getPostProbability(m));
                    likTmp4 *= contrib;
                  }
                totLikTmp4 += likTmp4;
              }
            }
          }
        }
        totLikTmp2 += likTmp2;
        totLikTmp3 += likTmp3;
        Message.debug("tot2="+totLikTmp2+" tot3="+totLikTmp3);
        totTotLikTmp4 += totLikTmp4;
      }
    }
    Message.debug("likTmp1="+likTmp1);
    Message.debug("totLikTmp2="+totLikTmp2);
    Message.debug("totLikTmp3="+totLikTmp3);
    Message.debug("totTotLikTmp4="+totTotLikTmp4);
    Message.debug("lik3a="+(likTmp1 + totLikTmp2 + totLikTmp3 + totTotLikTmp4));
    return likTmp1 + totLikTmp2 + totLikTmp3 + totTotLikTmp4;
  }

  /**
   * Computation for the carriers (of the ancestral haplotype): Group1
   * @param input       the input data
   * @param side        LEFT:0 RIGHT:1
   * @param nGeneration number of Generation
   * @param n          number of individuals belonging to the G1 group
   * @return the likelihood
   */
  private static double likeB(InputData input, int side, int nGeneration, int n) {
    int nMarker = input.getNMarker(side);
    double likB = 0;

    /* Loop on individual i, the ancestor */
    for (int i = 0; i < input.getNIndividuals(); i++) {
      final int endMarkerI = input.getEndMarker(side,i);
      //last marker of individual I is the last marker of the side
      if (endMarkerI == nMarker) {
        final double ui = C.U(input.getMutationRate(), endMarkerI, nGeneration);
        final double si = C.S(endMarkerI, nGeneration, input.getFractions(side));
        final double contribAnc = Math.pow(ui * si, input.getPostProbability(i));

        double likTmp1 = 1;
        double likTmp2 = 0;
        double likTmp3 = 0;
        double totLikTmp2 = 0;
        double totLikTmp3 = 0;
        double totTotLikTmp4 = 0;

        /* Loop on individual j */
        for (int j = 0; j < input.getNIndividuals(); j++) {
          double totLikTmp4 = 0.0;
          final int endMarkerJ = input.getEndMarker(side,j);
          //last marker of individual J is also the last marker of the side
          if ((endMarkerJ == nMarker) && (j != i)) {
            final double p1 = (endMarkerJ > 2) ? input.getFrequencies(side, endMarkerJ) : 0.0;
            final double p2 = (endMarkerJ > 2) ? input.getFrequencies(side, endMarkerJ - 1) : 0.0;
            final double u = C.U(input.getMutationRate(), endMarkerJ - 1, nGeneration);
            final double hj = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), input.getEndAlleles(side,i), input.getEndAlleles(side,j), input.isUseStepWiseModel());
            final double contrib = Math.pow(u * hj, input.getPostProbability(j));
            likTmp1 *= contrib;

            //If more than 1 G1
            if (n > 1) {

              final double f1 = C.F(endMarkerJ - 1, nGeneration, input.getFractions(side));
              final double f2 = C.F(endMarkerJ - 2, nGeneration, input.getFractions(side));
              //only last value will be kept
              likTmp2 = Math.pow(u * p1      * f1, input.getPostProbability(j));
              likTmp3 = Math.pow(u * p1 * p2 * f2, input.getPostProbability(j));
              for (int k = 0; k < input.getNIndividuals(); k++) {
                if ((input.getEndAlleles(side, k) == nMarker) && (k != i) && (k != j)) {
                  final double hk = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), input.getEndAlleles(side, i), input.getEndAlleles(side, k), input.isUseStepWiseModel());
                  final double tmpPow = Math.pow(u * hk, input.getPostProbability(k));
                  likTmp2 *= tmpPow;
                  likTmp3 *= tmpPow;
                }
              }
              //If more than 2 G1
              if (n > 2) {
                for (int l = j + 1; l < input.getNIndividuals(); l++) {
                  if ((input.getEndAlleles(side, l) == nMarker) && (l != i)) {
                    final double tmpCF = Math.pow(u * p1 * f1, 2);
                    double likTmp4 = Math.pow(tmpCF, input.getPostProbability(l) * input.getPostProbability(j));
                    for (int m = 0; m < input.getNIndividuals(); m++)
                      if ((input.getEndAlleles(side, m) == nMarker) && (m != l) && (m != i) && (m != j)) {
                        final double hm = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), input.getEndAlleles(side, i), input.getEndAlleles(side, m), input.isUseStepWiseModel());
                        final double tmpPow = Math.pow(u * hm, input.getPostProbability(m));
                        likTmp4 *= tmpPow;
                      }
                    totLikTmp4 += likTmp4;
                  }
                }
              }
            }

            totLikTmp2 += likTmp2;
            totLikTmp3 += likTmp3;
            totTotLikTmp4 += totLikTmp4;
          }
        }
        likB += contribAnc * (likTmp1 + totLikTmp2 + totLikTmp3 + totTotLikTmp4);
      }
    }
    return likB;
  }

  /**
   * Computation for the non carrier (of the ancestral haplotype) : Group 2
   * @param input       the input data
   * @param side        LEFT:0 RIGHT:1
   * @param nGeneration number of Generation
   * @return the likelihood
   */
  private static double likeC(InputData input, int side, int nGeneration) {
    int nMarker = input.getNMarker(side);
    /* Computations for individuals belonging to the G2 group */
    double likC = 1;
    for (int j = 0; j < input.getNIndividuals(); j++) {
      final int endMarkerJ = input.getEndMarker(side, j);
      //end marker of individual J is not the last marker of the side
      if (endMarkerJ != nMarker) {
        final double p1 = (endMarkerJ > 2) ? input.getFrequencies(side, endMarkerJ): 0.0;
        final double p2 = (endMarkerJ > 2) ? input.getFrequencies(side, endMarkerJ - 1) : 0.0;
        final double u = C.U(input.getMutationRate(), endMarkerJ - 1, nGeneration);
        final double h = C.H(input.getMutationRate(), endMarkerJ, nGeneration, input.getFractions(side), input.getAncestralAlleles(side,endMarkerJ), input.getEndAlleles(side,j), input.isUseStepWiseModel());
        final double f1 = C.F(endMarkerJ - 1, nGeneration, input.getFractions(side));
        final double f2 = C.F(endMarkerJ - 2, nGeneration, input.getFractions(side));
        final double contrib = u * h;
        likC *= Math.pow(contrib + p1 * f1 + p1 * p2 * f2, input.getPostProbability(j));
      }
    }
    return likC;
  }



  /**
   * Computation for subjects sharing the longest haplotype
   * @param input the input data
   * @param ng number of generation
   * @param nLongest number of individual in G1[left,right]
   * @return
   */
  private static double totLike(InputData input, int ng, final int[] nLongest) {
    double likeLeft  = input.hasLeft  ? like(input, C.LEFT,  ng, nLongest[C.LEFT])  : 1;
    double likeRight = input.hasRight ? like(input, C.RIGHT, ng, nLongest[C.RIGHT]) : 1;
    Message.debug("Left likelihood: "+likeLeft);
    Message.debug("Right likelihood: "+likeRight);
    return likeLeft * likeRight;
  }

  /**
   * Main computation
   * @param input the input data
   * @param nLongest number of individual in G1[left,right]
   * @return results of Estiage
   * @throws EstiageException if the maximum number of iterations is reached
   */
  private static Results findMaxLike(InputData input, final int[] nLongest) throws EstiageException {
    // Number of generations from common ancestor
    int nGen = 1;
    // Sum of the likelihood over nGen
    double pTot = 0;
    // Temporary variable to compute 95% CI
    double pTic = 0;
    // Keep the maximum likelihood at each incrementation of nGen
    double pMax = totLike(input, nGen, nLongest);
    // Array to keep the likelihood of n
    double[] pGenTot = new double[C.MAX_GENERATIONS + 1];
    int nMax = nGen;
    for (int gen = 1; gen <= C.MAX_GENERATIONS; gen++)
      pGenTot[gen] = 0.0;
    // Computation of the maximum likelihood<
    pGenTot[nGen] = totLike(input, nGen, nLongest);
    while (pMax / pGenTot[nGen] < C.LIKELIHOOD_STOP_RATIO) {
      if (pGenTot[nGen] >= pMax) {
        pMax = pGenTot[nGen];
        nMax = nGen;
      }
      nGen++;
      if (nGen > C.MAX_GENERATIONS) {
        String message = "Maximum number of iterations ["+nGen+"] reached";
        for (int i = 0; i < input.getNIndividuals(); i++)
          message += "\n" + i + " " + input.getEndMarker(C.LEFT,i) + " " + input.getEndMarker(C.RIGHT,i);
        throw new EstiageException(message);
      }
      pGenTot[nGen] = totLike(input, nGen, nLongest);
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

    return new Results(nMax, nEnd, nInf, nSup, lnTotLike);
  }
}

