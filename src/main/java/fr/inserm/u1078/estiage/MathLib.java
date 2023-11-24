package fr.inserm.u1078.estiage;

/**
 * Mathematical Functions
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-10-24
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class MathLib {
  //rate from hapmap is in cm/Mb, so cm = rate * dist_Mb
  //1cM = 1% recombination
  //tau is the recombination frequency in Morgans
  //theta is the observed recombination fraction
  //(1) tau = .25 * ln(1+2theta / 1-2theta)
  //(2) theta = (e^(4tau) - 1)/2(e^(4tau) + 1)

  /**
   * Compute tau - the frequency of genetic recombination in centimorgans between two positions
   * @param rRate the recombination rate between the positions
   * @param mBase the distance in Mb between the positions
   * @return
   */
  public static final double getCentiMorgan(double rRate, double mBase){
    return rRate * mBase;
  }

  /**
   * Uses the Kosambi formula to compute theta from tau
   * @param tau the frequency of genetic recombination in centimorgans
   * @return the recombination fraction theta
   */
  public static final double kosambiTheta(double tau){
    double exp = Math.exp(4*tau);
    if(Double.isInfinite(exp))
      return 0.5;
    double theta = 0.5 * (exp - 1) / (exp + 1);
    return theta;
  }

  /**
   * Uses the Kosambi formula to compute Tau from theta
   * @param theta the recombination fraction theta
   * @return tau - the frequency of genetic recombination in centimorgans
   */
  public static final double kosambiTau(double theta){
    //d=(1/4)ln((1+2r)/(1-2r))

    double num = 1+2*theta;
    double denom = 1-2*theta;
    return 0.25*Math.log(num / denom);
  }

  /**
   * Compute n!<br/>
   * Return type is double: can compute up to 170! after that returns infinity<br/>
   * with long: can compute up to 65! after that, returns negative or 0
   * @param n an integer
   * @return factorial(n)
   */
  public static final double fact(int n) {
    double prod = 1;
    //fact(0) = fact(1) = 1;
    for (int i = 2; i <= n; i++)
      prod *= i;
    return prod;
  }
}
