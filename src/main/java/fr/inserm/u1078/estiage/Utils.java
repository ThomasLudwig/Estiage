package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

/**
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-16
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Utils {

  /**
   * tabix command
   */
  public static final String DEFAULT_TABIX = "tabix";

  /**
   * Gets all the lines that cover a position from a tabixed VCF file
   * @param vcfFilename the name of the VCF file
   * @param chr the chromosome
   * @param position the position
   * @return All the matching lines from the VCF file
   * @throws IOException
   * @throws InterruptedException
   */
  public static ArrayList<String> getLinesFromTabixedVCF(String vcfFilename, String chr, int position) throws IOException, InterruptedException {
    return getLinesFromTabixedVCF(vcfFilename, chr, position, position);
  }

  /**
   * Gets all the lines that cover a position from a tabixed VCF file
   * @param vcfFilename the name of the VCF file
   * @param chr the chromosome
   * @param start the start position
   * @param end the end position
   * @return All the matching lines from the VCF file
   * @throws IOException
   * @throws InterruptedException
   */
  public static ArrayList<String> getLinesFromTabixedVCF(String vcfFilename, String chr, int start, int end) throws IOException, InterruptedException {
    return getLinesFromTabixedVCF(vcfFilename, chr+":"+start+"-"+end);
  }

  /**
   * Gets all the lines that cover a position from a tabixed VCF file
   * @param vcfFilename the name of the VCF file
   * @param pattern the query pattern ("chr", "chr:start-end")
   * @return All the matching lines from the VCF file
   * @throws IOException
   * @throws InterruptedException
   */
  public static ArrayList<String> getLinesFromTabixedVCF(String vcfFilename, String pattern) throws IOException, InterruptedException {
    ArrayList<String> ret = new ArrayList<>();

    String tabix = DEFAULT_TABIX;
    String overridetabix = System.getProperty("tabix");
    if(overridetabix != null && !overridetabix.isEmpty())
      tabix = overridetabix;
    String[] command = {tabix, vcfFilename, pattern};

    ProcessBuilder pb = new ProcessBuilder(command);
    pb.redirectErrorStream(true);
    Process process = pb.start();
    BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
    String line;
    while((line = in.readLine()) != null)
      ret.add(line);
    BufferedReader err = new BufferedReader(new InputStreamReader(process.getErrorStream()));
    String errline;
    String message = "";
    while((errline = err.readLine()) != null)
      message += errline + "\t";
    process.waitFor();

    in.close();
    err.close();
    if(!message.isEmpty())
      Message.error(message);
    Message.info("Found "+ret.size());
    return ret;
  }
}
