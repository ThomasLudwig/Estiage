package fr.inserm.u1078.estiage.ctranslation;

/**
 * Class representing the Results of Estiage.
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-20
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Results {
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
