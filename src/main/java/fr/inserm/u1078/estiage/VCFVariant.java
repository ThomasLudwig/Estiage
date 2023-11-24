package fr.inserm.u1078.estiage;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * A Variant within a VCF file
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-13
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class VCFVariant {
  private final String chr;
  private final int pos;
  private final String id;
  private final String[] alleles;
  private final String[] genotypes;

  public static final String MISSING = "-1";
  public static final String HETERO = "-2";

  //TODO rewrite entirely to work on phased data or write a wrapper that produce a HOMO file from phased data

  /**
   * Creates a Variant from a VCF line
   * @param line
   */
  public VCFVariant(String line) {
    String[] f = line.split("\t");
    this.chr = f[0];
    this.pos = Integer.parseInt(f[1]);
    this.id = (f[2].isEmpty() || ".".equals(f[2])) ? chr+":"+pos : f[2];
    this.alleles = (f[3]+","+f[4]).split(",");
    this.genotypes = new String[f.length - 9];
    for(int i = 0; i < this.genotypes.length; i++)
      this.genotypes[i] = getHomozygousAlleleFromGenotype(f[i+9]);
  }

  /**
   * check if a given variant is in the current VCF file
   * @param chr the chromosome to consider
   * @param pos the position to consider
   * @param allele the allele to consider
   * @param mode the search mode
   * @return true if the given variant is present for the given search mode
   */
  public boolean has(String chr, int pos, String allele, VCFFile.Mode mode){
    if(!this.chr.equalsIgnoreCase(chr))
      return false;
    if(this.pos != pos)
      return false;
    if(mode == VCFFile.Mode.IGNORE)
      return true;
    for(String al : alleles) {
      if (al.equalsIgnoreCase(allele) && mode == VCFFile.Mode.HOMOZYGOUS)
        return true;
      if (al.equalsIgnoreCase(HETERO) && mode == VCFFile.Mode.HETEROZYGOUS)
        return true;
    }
    return false;
  }

  /**
   * Gets the genotype as a String
   * @param s
   * @return the homozygous allele (0,1,2...), -1 if it as a missing genotype, -2 if the genotype is heterozygous
   */
  private String getHomozygousAlleleFromGenotype(String s) {
    if(s.startsWith("."))
      return MISSING;
    String[] f = s.split(":");
    String geno = f[0].replace("|","/");
    String[] g = geno.split("/");
    int g1 = Integer.parseInt(g[0]);
    int g2 = Integer.parseInt(g[1]);
    if(g1 != g2)
      return HETERO;
    return alleles[g1];
  }

  /**
   * checks if the given allele can be a target in the given mode
   * @param targetAllele
   * @param mode
   * @return
   */
  public boolean canBeTarget(String targetAllele, VCFFile.Mode mode){
    for(String geno : genotypes) {
      if (!targetAllele.equalsIgnoreCase(geno) && mode == VCFFile.Mode.HOMOZYGOUS)
        return false;//all allele should be the same
      if (!targetAllele.equalsIgnoreCase(HETERO) && mode == VCFFile.Mode.HETEROZYGOUS)
        return false;//all allele should be hetero
    }
    return true;
  }

  /**
   * Checks if the Variant is valid
   * @return false if at least one genotype is HETERO or Missing
   */
  public boolean isValid(){
    for(String geno : genotypes){
      if(HETERO.equals(geno))
        return false;//allele should be homozygous
      if(MISSING.equals(geno))
        return false;//no missing allele
    }
    return true;
  }

  /**
   * Gets the genotypes for the variant
   * @return
   */
  public String[] getGenotypes(){
    return genotypes;
  }

  /**
   * Gets the chromosome
   * @return
   */
  public String getChr() {
    return chr;
  }

  /**
   * Gets the position of the variant
   * @return
   */
  public int getPos() {
    return pos;
  }

  /**
   * Gets the ID of the variant
   * @return
   */
  public String getId() {
    return id;
  }

  /**
   * Gets the REF/ALT alleles
   * @return
   */
  public String[] getAlleles() {
    return alleles;
  }

  /**
   * TopAlleles and the number of samples that have this allele
   * Examples : {"14", "A"} or {"4","A","G"} or {"1","A","C","G"}
   * @param samples the list of sample to take into account
   * @return a String[] containing the list of TopAlleles and the count (at the 0 index).
   */
  public String[] getTopAllelesAndCount(ArrayList<Integer> samples) {
    //Get the count for each allele
    HashMap<String, Integer> counts = new HashMap<>();
    for(int i : samples)
      counts.merge(genotypes[i], 1, Integer::sum);

    //Sort alleles by count DESC
    //number of allele
    //since there will not be a lot of values, this code doesn't have to be spotless
    int n = counts.size();
    String[] alleles = new String[n];
    int[] occurrences = new int[n];
    //n times, fill alleles/occurrences arrays and remove the allele
    for(int i = 0 ; i < n; i++){
      int max = -1;
      String maxAllele="-";
      for(String allele : counts.keySet()){
        int count = counts.get(allele);
        if(count > max){
          max = count;
          maxAllele = allele;
        }
      }
      alleles[i] = maxAllele;
      occurrences[i] = max;
      counts.remove(maxAllele);
    }
    int max = occurrences[0];
    String ret = max+","+alleles[0];
    for(int i = 1 ; i < n; i++) {
      if (occurrences[i] == max)
        ret += "," + alleles[i];
      else
        break;
    }
    return ret.split(",");
  }

  /**
   * Removes the samples that do not carry one of the TOP alleles (unique top allele = ancestral allele)
   * @param samples the list of samples to test
   * @param topAlleles the description of top alleles NumberOfSamples:Allele1:Allele2:.../AlleleN
   */
  public void removeNonAncestral(ArrayList<Integer> samples, String[] topAlleles){
    //If more there one top allele, clear for everyone
    if(topAlleles.length > 2) {
      samples.clear();
      return;
    }

    //else, clear only samples with another allele
    String ancestral = topAlleles[1] ;
    ArrayList<Integer> toRemove = new ArrayList<>();
    for(int i : samples)
      if(!ancestral.equals(genotypes[i]))
        toRemove.add(i);
    samples.removeAll(toRemove);
  }
}
