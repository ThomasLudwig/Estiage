package fr.inserm.u1078.estiage.test;

import fr.inserm.u1078.estiage.ctranslation.Estiage;
import fr.inserm.u1078.estiage.MathLib;
import fr.inserm.u1078.estiage.TSVFile;

/**
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-05-03
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class MainTest {
  public static void main(String[] arg) throws Exception {
    //phase();
    //testKosambi();

    int N = 200;
    long lfact = 1;
    double dfact = 1;
    for(int i = 2; i <= N; i++){
      lfact *= i;
      dfact *= i;
      System.out.println(i+" -> "+lfact+" / "+dfact);
    }
  }

  private static void testKosambi(){
    double tau = 0.2776; //in Morgan
    for(int i = 1 ; i <= 10; i++){
      double theta = MathLib.kosambiTheta(tau);
      double tau2 = MathLib.kosambiTau(theta);
      System.out.println("Tau "+tau+"->"+tau2+"     "+(tau/tau2)+" ----------- "+theta);
    }
  }

  private static final void phase() throws Exception {
    String input = "C:\\Users\\Thomas Ludwig\\Documents\\Projet\\Estiage\\data.unphased.teresa.txt";
    fr.inserm.u1078.estiage.Main.main(new String[]{fr.inserm.u1078.estiage.Main.KEY_PHASE, input, "10", "chr19:4567894"});
  }

  private static final void testLoad() throws Exception {
    String input = "C:\\Users\\user\\Downloads\\raw20220503.tsv";

    TSVFile tsv = new TSVFile(input, TSVFile.Type.RAW);
    System.out.println("Samples : "+tsv.getSamples().length);
    System.out.println(String.join("\t", tsv.getSamples()));
    System.out.println("Left : "+tsv.getLeftMarkers().length);
    System.out.println("Right : "+tsv.getRightMarkers().length);
  }

  private static final void testrun() throws Exception{
    String dat = "C:\\Users\\user\\Downloads\\emmanuelle.mathilde.dat";
    Estiage.run(dat);
  }
}
