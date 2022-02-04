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

  public GnomAD(String filename) throws IOException {
    this.filename = filename;
    this.tabixFilename = filename + ".tbi";
    this.check();
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
    this.hasChr = hasChr;
  }

  private void check() throws IOException {
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
  }

  public void applyFrequency(Marker m) throws IOException, EstiageFormatException {
    Message.info("Looking in ["+this.filename+"] for ["+m.getChromosome()+":"+m.getPosition()+"]");
    try {
      double frq = getFrequency(m.getChromosome(), m.getPosition(), m.getAncestral());
      m.setFrequencies(frq);
    } catch(InterruptedException e) {
      Message.error("InterrupedException ["+e.getMessage()+"] while looking for frequency of marker ["+m.getName()+"]");
    }
  }

  public double getFrequency(String chr, int position, String allele) throws IOException, InterruptedException, EstiageFormatException {
    for(String line : getLines(chr, position)) {
      String[] f = line.split("\t", -1);
      //If position was found //TODO what of ACT->ACG in N-2 ?
      if(position == Integer.parseInt(f[1])){
        //Search for the correct alt
        int idx = -1;
        String[] alleles = f[4].split(",");
        for(int i = 0 ; i < alleles.length; i++)
          if(alleles[i].equals(allele)) {
            idx = i;
            break;
          }

        //if alt was found, get its AF
        if(idx > -1) {
          String[] info = f[7].split(";", -1);
          for (String inf : info) {
            if (inf.startsWith("AF=")) {
              String[] afs = inf.substring(3).split(",", -1);
              try{
                return Double.parseDouble(afs[idx]);
              } catch(NumberFormatException e){
                return 0;
              }
            }
          }
        }
      }
    }
    return 0;
  }

  public String findChromosomes(String chr) throws EstiageFormatException {
    for(String[] chrs : CHROMOSOMES){
      if(chrs[0].equalsIgnoreCase(chr) || chrs[1].equalsIgnoreCase(chr)){
        int idx = hasChr ? 1 : 0;
        return chrs[idx];
      }
    }
    throw new EstiageFormatException("Unknown chromosome ["+chr+"]");
  }

  public ArrayList<String> getLines(String chr, int position) throws IOException, InterruptedException, EstiageFormatException {
    return subGetLines(findChromosomes(chr), position);
  }

  public ArrayList<String> subGetLines(String chr, int position) throws IOException, InterruptedException {
    ArrayList<String> ret = new ArrayList<>();


    //tabix this.filename chr:position-position > outvcf;
    String[] command = {"/PROGS/bin/tabix", filename, chr+":"+position+"-"+position};

    //NEW
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
