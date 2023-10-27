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

  public static final double getCentiMorgan(double rRate, double mBase){
    return rRate * mBase;
  }

  public static final double kosambiTheta(double tau){
    //rate from hapmap is in cm/Mb, so cm = rate * dist_Mb
    //1cM = 1% recombination
    //d is the recombination frequency in Morgans
    //r is the observed recombination fraction
    //(1) d = .25 * ln(1+2r / 1-2r)
    //(2) r = (e^(4d) - 1)/2(e^(4d) + 1)
    double exp = Math.exp(4*tau);
    if(Double.isInfinite(exp))
      return 0.5;
    double theta = 0.5 * (exp - 1) / (exp + 1);
    return theta;
  }

  public static final double kosambiTau(double theta){
    //d=(1/4)ln((1+2r)/(1-2r))

    double num = 1+2*theta;
    double denom = 1-2*theta;
    return 0.25*Math.log(num / denom);
  }

  /**
   * Compute n!
   * @param n an integer
   * @return factorial(n)
   */
  public static final double fact(int n) {
    double prod = 1.0;
    //fact(0) = fact(1) = 1;
    for (int i = 2; i <= n; i++)
      prod *= i;
    return prod;
  }
}
