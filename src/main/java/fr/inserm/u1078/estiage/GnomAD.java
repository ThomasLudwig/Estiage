package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.*;
import java.util.ArrayList;

/**
 * Class used to request GnomAD frequencies
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-01-11
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class GnomAD {
  private final String filename;
  private final String tabixFilename;
  private final boolean hasChr;

  //The name of the chromosomes without/with chr. Can't use numbers because of tabix
  public static final String[][] CHROMOSOMES = {
          {"1", "chr1"},
          {"2", "chr2"},
          {"3", "chr3"},
          {"4", "chr4"},
          {"5", "chr5"},
          {"6", "chr6"},
          {"7", "chr7"},
          {"8", "chr8"},
          {"9", "chr9"},
          {"10", "chr10"},
          {"11", "chr11"},
          {"12", "chr12"},
          {"13", "chr13"},
          {"14", "chr14"},
          {"15", "chr15"},
          {"16", "chr16"},
          {"17", "chr17"},
          {"18", "chr18"},
          {"19", "chr19"},
          {"20", "chr20"},
          {"21", "chr21"},
          {"22", "chr22"},
          {"X", "chrX"},
          {"Y", "chrY"},
          {"MT", "chrM"}
  };

  /**
   * Creates a new GnomAD object and checks if the provided file is valid
   * @param filename the name of the GnomAD file
   * @throws IOException
   */
  public GnomAD(String filename) throws IOException {
    this.filename = filename;
    this.tabixFilename = filename + ".tbi";
    this.hasChr = this.check();
  }

  /**
   * Checks that <ul>
   *   <li>the vcf exists</li>
   *   <li>the vf is not a directory</li>
   *   <li>the vcf is bgzipped</li>
   *   <li>the tabix exists</li>
   *   <li>the tabix is not a directory</li>
   * </ul>
   * @return true if the chromosome names start with "chr"
   * @throws IOException
   */
  private boolean check() throws IOException {
    File vcf = new File(filename);
    File tabix = new File(tabixFilename);
    if(!vcf.exists())
      throw new FileNotFoundException("File "+filename+" does not exist");
    if(vcf.isDirectory())
      throw new FileNotFoundException("File "+filename+" is a directory");
    if(!filename.toLowerCase().endsWith(".gz"))
      throw new IOException("File "+filename+" does not seem to be bgzipped");
    if(!tabix.exists())
      throw new FileNotFoundException("File "+filename+" does not exist");
    if(tabix.isDirectory())
      throw new FileNotFoundException("File "+filename+" is a directory");

    UniversalReader in = new UniversalReader(filename);
    String line;
    boolean hasChr = false;
    while((line = in.readLine()) != null){
      if(!line.startsWith("#")){
        hasChr = line.startsWith("chr");
        break;
      }
    }
    in.close();
    return hasChr;
  }

  /**
   * Get the allele frequency from the GnomAD file
   * @param chr the chromosome
   * @param position the position
   * @param allele the allele
   * @return the allele frequency
   * @throws IOException
   * @throws EstiageFormatException
   */
  public double getFrequency(String chr, int position, String allele) throws IOException, EstiageFormatException {
    try {
      String tabixChr = findChromosomes(chr);
      for (String line : Utils.getLinesFromTabixedVCF(this.filename, tabixChr, position)) {
        String[] f = line.split("\t", -1);
        //If position was found //TODO what of ACT->ACG in N-2 ?
        if (position == Integer.parseInt(f[1])) {
          //Search for the correct alt
          int idx = -1;
          String[] alleles = f[4].split(",");
          for (int i = 0; i < alleles.length; i++)
            if (alleles[i].equals(allele)) {
              idx = i;
              break;
            }

          //if alt was found, get its AF
          if (idx > -1) {
            String[] info = f[7].split(";", -1);
            for (String inf : info) {
              if (inf.startsWith("AF=")) {
                String[] afs = inf.substring(3).split(",", -1);
                try {
                  return Double.parseDouble(afs[idx]);
                } catch (NumberFormatException e) {
                  return 0;
                }
              }
            }
          }
        }
      }
    }
    catch(InterruptedException e){
      Message.error("InterrupedException ["+e.getMessage()+"] while looking for frequency for ["+chr+":"+position+":"+allele+"]");
    }
    return 0;
  }

  /**
   * Return the chromosome name known in the tabix file for a given chromosome
   * @param chr the choromose name in hg/GRCh format
   * @return
   * @throws EstiageFormatException
   */
  public String findChromosomes(String chr) throws EstiageFormatException {
    for(String[] chrs : CHROMOSOMES){
      if(chrs[0].equalsIgnoreCase(chr) || chrs[1].equalsIgnoreCase(chr)){
        int idx = hasChr ? 1 : 0;
        return chrs[idx];
      }
    }
    throw new EstiageFormatException("Unknown chromosome ["+chr+"]");
  }
}
