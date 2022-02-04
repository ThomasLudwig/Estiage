package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.util.Collection;
import java.util.HashMap;

/**
 * Representation of a Marker within a TSVFile
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2021-03-18
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Marker {
  public static final String MISSING = "";

  private final String name;
  private String chr;
  private int position;
  private double megaBases;
  private double distance;
  private double cM;
  private double recombinationFraction;
  private double rate;
  private double frequencies;
  private final HashMap<String, String> haplotypes;
  private String ancestral = "";

  public Marker(String name) {
    this.name = name;
    this.haplotypes = new HashMap<>();
  }

  public Marker(VCFFile.Variant variant, String[] samples){
    this.name = variant.getId();
    this.chr = variant.getChr();
    this.position = variant.getPos();
    haplotypes = new HashMap<>();
    for(int i = 0 ; i < samples.length; i++)
      haplotypes.put(samples[i], variant.getGenotypes()[i]);
  }

  public Marker(String chrposallele, String[] samples){
    String[] f = chrposallele.split(":");
    this.name = chrposallele;
    this.chr = f[0];
    this.position = Integer.parseInt(f[1]);
    haplotypes = new HashMap<>();
    for(int i = 0 ; i < samples.length; i++)
      haplotypes.put(samples[i], f[2]);
  }

  public void setDistance(Marker target){
    this.distance = Math.abs(this.megaBases - target.megaBases);
  }

  public void setChromosomeAndPosition(String chrpos){
    try{
      String[] f = chrpos.split(":");
      this.setChromosome(f[0]);
      this.setPosition(Integer.parseInt(f[1]));
    } catch (ArrayIndexOutOfBoundsException | NumberFormatException e){
      Message.error("Unexpected format for ["+chrpos+"], should be chr:pos");
    }
  }

  public void setPosition(int position){
    this.position = position;
    this.megaBases = this.position / 1000000D;
  }

  public void setChromosome(String chr){
    this.chr = chr;
  }

  public void addHaplotype(String sample, String haplotype){
    this.haplotypes.put(sample, haplotype);
  }

  /**
   *
   */
  public void setAncestral(){
    //Counts frequency
    HashMap<String, Integer> counts = new HashMap<>();
    for(String val : haplotypes.values()){
      if(val != MISSING) {
        Integer nb = counts.get(val);
        if (nb == null)
          nb = 0;
        counts.put(val, nb + 1);
      }
    }
    //get most frequent counts
    int max = 0;
    for(Integer count : counts.values())
      if(count > max)
        max = count;
    //get most frequent alleles, if unique return
    String allele = "-1";
    for(String a : counts.keySet()){
      if(counts.get(a) == max){
        if("-1".equals(allele)) {
          allele = a;
        } else {
          //if ex aequo
          allele = "-1";
          break;
        }
      }
    }
    this.ancestral = allele;
  }

  public String getAncestral(){
    return this.ancestral;
  }

  public boolean isEmpty(){
    for(String val : haplotypes.values())
      if(val != null)
        return false;
    return true;
  }

  public int getPosition(){
    return this.position;
  }

  public void setFrequencies(double frequencies) {
    this.frequencies = frequencies;
  }

  public void setRate(double rate){
    this.rate = rate;
    applyKosambi();
  }

  /**
   * kosambi function
   * See http://www.crypticlineage.net/lib/Kosambi.pdf
   * https://www.ias.ac.in/public/Volumes/reso/016/06/0540-0550.pdf page 9
   */
  public void applyKosambi() {
    //rate from hapmap is in cm/Mb, so cm = rate * dist_Mb
    //1cM = 1% recombination
    this.cM = this.rate * this.distance;
    //d is the recombination frequency in Morgans
    //r is the observed recombination fraction
    double d = this.cM / 100; //from cMorgans to Morgans
    //(1) d = .25 * ln(1+2r / 1-2r)
    //(2) r = (e^(4d) - 1)/2(e^(4d) + 1)
    double exp = Math.exp(4*d);
    this.recombinationFraction = 0.5 * (exp - 1) / (exp + 1);
  }

  public String getChromosome(){
    return chr;
  }

  public String getName() {
    return name;
  }

  public double getDistance() {
    return distance;
  }

  public double getMegaBases() {
    return megaBases;
  }

  public void setMegaBases(double megaBases) {
    this.megaBases = megaBases;
  }

  public double getcM() {
    return cM;
  }

  public double getRecombinationFraction() {
    return recombinationFraction;
  }

  public void setRecombinationFraction(double recombinationFraction) { this.recombinationFraction = recombinationFraction; }

  public double getRate() {
    return rate;
  }

  public double getFrequencies() { return frequencies; }

  public String getHaplotype(String sample) {
    return haplotypes.get(sample);
  }

  public void setHaplotype(String sample, String haplotype) {
    this.haplotypes.put(sample, haplotype);
  }

  public Collection<String> getAllAlleles(){
    return this.haplotypes.values();
  }
}