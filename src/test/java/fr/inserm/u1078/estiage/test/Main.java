package fr.inserm.u1078.estiage.test;

import fr.inserm.u1078.estiage.EstiageCTranslation;
import fr.inserm.u1078.estiage.TSVFile;

/**
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-05-03
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Main {
  public static void main(String[] arg) throws Exception {
    testrun();
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
    EstiageCTranslation e = new EstiageCTranslation(dat);
    e.run();
  }
}
