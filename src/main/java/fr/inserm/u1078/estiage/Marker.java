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
  private double distanceMb;
  private double cM;
  private double recombinationFraction;
  private double rate;
  private double frequency;
  private final HashMap<String, String> alleles;
  private String ancestral = "";

  /**
   * Creates an empty Marker with the given name
   * @param name the name of the marker
   */
  public Marker(String name) {
    this.name = name;
    this.alleles = new HashMap<>();
  }

  /**
   * Creates a Marker for the given variant and samples
   * @param variant the variant
   * @param samples the samples
   */
  public Marker(VCFVariant variant, String[] samples){
    this.name = variant.getId();
    this.chr = variant.getChr();
    this.position = variant.getPos();
    alleles = new HashMap<>();
    for(int i = 0 ; i < samples.length; i++)
      alleles.put(samples[i], variant.getGenotypes()[i]);
  }

  /**
   * Creates a Marker for the given allele and samples
   * @param chrposallele the variant in the format chromosome:position:allele
   * @param samples the samples
   */
  public Marker(String chrposallele, String[] samples){
    String[] f = chrposallele.split(":");
    this.name = chrposallele;
    this.chr = f[0];
    this.position = Integer.parseInt(f[1]);
    alleles = new HashMap<>();
    for(int i = 0 ; i < samples.length; i++)
      alleles.put(samples[i], f[2]);
  }

  /**
   * Set the distance in Megabase between this marker and the target marker
   * @param target the target marker
   */
  public void setDistanceMb(Marker target){
    this.distanceMb = Math.abs(this.megaBases - target.megaBases);
  }

  /**
   * sets the chromosome and position for this marker
   * @param chrpos chromosome:position
   */
  public void setChromosomeAndPosition(String chrpos){
    try{
      String[] f = chrpos.split(":");
      this.setChromosome(f[0]);
      this.setPosition(Integer.parseInt(f[1]));
    } catch (ArrayIndexOutOfBoundsException | NumberFormatException e){
      Message.error("Unexpected format for ["+chrpos+"], should be chr:pos");
    }
  }

  /**
   * sets the position for this marker
   * @param position
   */
  public void setPosition(int position){
    this.position = position;
    this.megaBases = this.position / 1000000D;
  }

  /**
   * sets the chromosome of this marker
   * @param chr
   */
  public void setChromosome(String chr){
    this.chr = chr;
  }

  /**
   * Compute the ancestral allele. If None is found, ancestral allele is set to the default missing value ""
   */
  public void setAncestral(){
    //Counts frequency
    HashMap<String, Integer> counts = new HashMap<>();
    for(String val : alleles.values()){
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

  /**
   * Gets the ancestral allele for this marker
   * @return
   */
  public String getAncestral(){
    return this.ancestral;
  }

  /**
   * checks of this marker is empty
   * @return
   */
  public boolean isEmpty(){
    for(String val : alleles.values())
      if(val != null)
        return false;
    return true;
  }

  /**
   * Gets the position of this marker
   * @return
   */
  public int getPosition(){
    return this.position;
  }

  /**
   * Sets the frequency for the ancestral allele of this marker
   * @param frequency the frequency
   */
  public void setFrequency(double frequency) {
    this.frequency = frequency;
  }

  /**
   * Sets the following values, between this marker and the target<ul>
   *   <li>recombination rate</li>
   *   <li>distance in cM</li>
   *   <li>recombination fraction</li>
   * </ul>
   * @param rate the recombination rate between this marker and the target
   */
  public void setRate(double rate){
    this.rate = rate;
    //rate from hapmap is in cm/Mb, so cm = rate * dist_Mb
    //1cM = 1% recombination
    this.cM = MathLib.getCentiMorgan(this.rate, this.distanceMb);
    double d = this.cM / 100; //from cMorgans to Morgans
    this.recombinationFraction = MathLib.kosambiTheta(d);
  }

  /**
   * gets the chromosome for this marker
   * @return
   */
  public String getChromosome(){
    return chr;
  }

  /**
   * gets the name of this marker
   * @return
   */
  public String getName() {
    return name;
  }

  /**
   * Gets the distance in Mb between this marker and the target
   * @return
   */
  public double getDistanceMb() {
    return distanceMb;
  }

  /**
   * Gets the position of this marker in Mb
   * @return
   */
  public double getMegaBases() {
    return megaBases;
  }

  /**
   * Sets the position of this marker in Mb
   * @param megaBases the position in Mb
   */
  public void setMegaBases(double megaBases) {
    this.megaBases = megaBases;
  }

  /**
   * Gets get distance in cM between this marker and the target
   * @return
   */
  public double getcM() {
    return cM;
  }

  /**
   * Gets the recombination fraction theta between this marker and the target
   * @return
   */
  public double getRecombinationFraction() {
    return recombinationFraction;
  }

  /**
   * sets the recombination fraction theta between this marker and the target
   * @param recombinationFraction theta
   */
  public void setRecombinationFraction(double recombinationFraction) { this.recombinationFraction = recombinationFraction; }

  /**
   * Gets the recombination rate between this marker and the target
   * @return
   */
  public double getRate() {
    return rate;
  }

  /**
   * Gets the frequency of the ancestral allele for this marker
   * @return
   */
  public double getFrequency() { return frequency; }

  /**
   * Gets the allele for the given sample
   * @param sample the name of the sample
   * @return
   */
  public String getAllele(String sample) {
    return alleles.get(sample);
  }

  /**
   * Sets the allele for the given sample
   * @param sample the sample name
   * @param allele the allele
   */
  public void setAllele(String sample, String allele) {
    this.alleles.put(sample, allele);
  }

  /**
   * Gets allele alleles for this marker
   * @return unsorted collection of allele
   */
  public Collection<String> getAllAlleles(){
    return this.alleles.values();
  }
}