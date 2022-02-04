package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

/**
 * VCF Parser, used to transform VCF Data into a TSV Raw file
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-01-20
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class VCFFile {

  public enum Mode {IGNORE, HETEROZYGOUS, HOMOZYGOUS}

  private final String filename;
  private final boolean isTabix;
  private final String[] samples;
  private final Mode mode;

  private String chr;
  private int pos;
  private String allele;

  public VCFFile(String filename, Mode mode) throws IOException, EstiageFormatException {
    this.filename = filename;
    this.mode = mode;
    this.isTabix = isTabix();
    samples = readSamples();
  }

  private boolean isTabix() throws FileNotFoundException {
    File vcf = new File(filename);
    File tabix = new File(filename+".tbi");
    if(!vcf.exists())
      throw new FileNotFoundException("File "+filename+" does not exist");
    if(vcf.isDirectory())
      throw new FileNotFoundException("File "+filename+" is a directory");
    if(!filename.toLowerCase().endsWith(".gz"))
      return false;
    if(!tabix.exists())
      return false;
    if(tabix.isDirectory())
      return false;
    return true;
  }

  private String[] readSamples() throws IOException, EstiageFormatException {
    UniversalReader in = new UniversalReader(this.filename);
    String line;
    String header = null;
    while((line = in.readLine()) != null){
      if(line.startsWith("#")) {
        if(line.startsWith("#CHROM"))
          header = line;
      }
      else
        break;
    }
    in.close();
    if(header == null)
      throw new EstiageFormatException("No header found in VCF ["+filename+"]");
    String[] f = header.split("\t");
    String[] s = new String[f.length - 9];
    if (s.length >= 0) System.arraycopy(f, 9, s, 0, s.length);
    return s;
  }

  public void setVariant(String chrPosAllele) {
    String[] f = chrPosAllele.split(":");
    this.chr = f[0];
    this.pos = Integer.parseInt(f[1]);
    this.allele = f[2];
  }

  private int[] getLeftRight(ArrayList<Variant> variants){
    int size = variants.size();
    int left = Integer.MAX_VALUE;
    int right = -1;

    int start = 0;
    int end = variants.size() - 1;

    while(start <= end ) {
      int i = (start + end) / 2;
      //Message.info("["+start+"|"+i+"|"+end+"]");
      int v = variants.get(i).pos;
      if (v == pos) {
        return new int[]{i - 1, i + 1};
      } else if (v < pos) {
        if (i == size - 1)
          return new int[]{i, i + 1};
        if (variants.get(i + 1).pos > pos)
          return new int[]{i, i + 1};
        start = i;
      } else { //v > pos
        if (i == 0)
          return new int[]{i - 1, i};
        if (variants.get(i - 1).pos < pos)
          return new int[]{i - 1, i};
        end = i;
      }
    }

    return new int[]{left, right};
  }

  public void exportAsRaw(String raw) throws IOException, EstiageFormatException, InterruptedException {
    ArrayList<Variant> variants = loadFile();
    Message.info(variants.size()+ " valid variants found in ["+filename+"]");
    int[] leftRight = getLeftRight(variants);
    int left = leftRight[0];
    int right = leftRight[1];
    //Message.info("Looking for ["+pos+"] and found here ["+left+","+right+"] ["+variants.get(left).pos+";"+variants.get(right).pos+"]");
    if(left > right)
      throw new EstiageFormatException("Something went wrong. Is you VCF file unsorted?");
    Marker target;
    if(left + 2 == right) {
      Variant targetVariant = variants.get(left + 1);
      //Message.info("Target is "+targetVariant.pos);
      if(!targetVariant.has(chr, pos, allele, mode))
        throw new EstiageFormatException("The genotypes for the variant at ["+chr+":"+pos+":"+allele+"] do not match the selected mode ["+mode+"]");
       target = new Marker(targetVariant, samples);
    } else {
      target = new Marker(chr+":"+pos+":"+allele, samples);
      //Message.info("Target is default");
    }
    /*
      int idx = -1;
      for(int i = 0 ; i < variants.size(); i++){
        if(variants.get(i).has(chr, pos, allele)){
          idx = i;
          break;
        }
      }
      if(idx == -1)
        throw new EstiageFormatException("Variant ["+chr+":"+pos+":"+allele+"] has not been found in the file ["+filename+"]");
     */
    Message.info("Target Variant is between  ["+(left+1)+";"+(right+1)+"]/"+variants.size());
    ArrayList<Integer> remainingSamples = new ArrayList<>();

    //Left
    TreeMap<Integer, ArrayList<Variant>> leftVariants = new TreeMap<>();
    for(int s = 0 ; s < samples.length; s++)
      remainingSamples.add(s);

    for(int i = left; i >= 0 && !remainingSamples.isEmpty(); i--){
      Variant v = variants.get(i);
      v.removeNonAncestral(remainingSamples);
      int size = remainingSamples.size();
      ArrayList<Variant> tmpVariants = leftVariants.get(size);
      if(tmpVariants == null)
        tmpVariants = new ArrayList<>();
      tmpVariants.add(v);
      leftVariants.put(size, tmpVariants);
    }

    //Right
    TreeMap<Integer, ArrayList<Variant>> rightVariants = new TreeMap<>();
    for(int s = 0 ; s < samples.length; s++)
      remainingSamples.add(s);
    for(int i = right ; i < variants.size() && !remainingSamples.isEmpty(); i++){
      Variant v = variants.get(i);
      v.removeNonAncestral(remainingSamples);
      int size = remainingSamples.size();
      ArrayList<Variant> tmpVariants = rightVariants.get(size);
      if(tmpVariants == null)
        tmpVariants = new ArrayList<>();
      tmpVariants.add(v);
      rightVariants.put(size, tmpVariants);
    }

    ArrayList<Marker> leftMarkers = new ArrayList<>();
    for(int smpls : leftVariants.descendingKeySet()){
      ArrayList<Variant> vars = leftVariants.get(smpls);
      leftMarkers.add(new Marker(vars.get(0), samples));
      if(vars.size() > 1)
        leftMarkers.add(new Marker(vars.get(vars.size() - 1), samples));
    }

    ArrayList<Marker> rightMarkers = new ArrayList<>();
    for(int smpls : rightVariants.descendingKeySet()){
      ArrayList<Variant> vars = rightVariants.get(smpls);
      rightMarkers.add(new Marker(vars.get(0), samples));
      if(vars.size() > 1)
        rightMarkers.add(new Marker(vars.get(vars.size() - 1), samples));
    }
    Message.info("Variants kept on the left side : "+leftMarkers.size());
    Message.info("Variants kept on the right side : "+rightMarkers.size());
    TSVFile tsvFile = new TSVFile(target, leftMarkers.toArray(new Marker[0]), rightMarkers.toArray(new Marker[0]), samples);
    tsvFile.export(raw);
  }

  private ArrayList<Variant> getVariantsFromTabix(String pattern, boolean onlyValid) throws IOException, InterruptedException {
    ArrayList<Variant> variants = new ArrayList<>();
    String[] command = {"/PROGS/bin/tabix", filename, pattern};

    ProcessBuilder pb = new ProcessBuilder(command);
    pb.redirectErrorStream(true);
    Process process = pb.start();
    BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
    String line;
    while((line = in.readLine()) != null) {
      Variant v = new Variant(line);
      if(!onlyValid || v.isValid())
        variants.add(v);
    }
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
    return variants;
  }

  private void checkHasVariant() throws EstiageFormatException, IOException, InterruptedException {
    if(this.mode == Mode.IGNORE)
      return;
    if(!this.isTabix) {
      Message.warning("As the fill is not tabixed, no check will be performed on the target variant's genotypes");
      return;
    }

    for(Variant v : this.getVariantsFromTabix(chr+":"+pos+"-"+pos, false)) {
      if (v.has(chr, pos, allele, mode)) {
        if (!v.canBeTarget(allele, mode))
          throw new EstiageFormatException("Variant to estimate should be homozygous for each sample in the file. Here: " + String.join(",", v.getGenotypes()));
        return; //Found
      }
    }
    throw new EstiageFormatException("VCF file [" + filename + "] doesn't not contains a variant for " + chr + ":" + pos + ":" + allele);
  }

  private ArrayList<Variant> loadFile() throws IOException, EstiageFormatException, InterruptedException {
    this.checkHasVariant();
    if(this.isTabix){
      return this.getVariantsFromTabix(chr, true);
    } else {
      ArrayList<Variant> variants = new ArrayList<>();
      Message.warning("File ["+filename+"] is not tabixed, this will be slow");
      UniversalReader in = new UniversalReader(this.filename);
      String line;
      int read = 0;
      while((line = in.readLine()) != null){
        if(line.startsWith(chr)){
          read++;
          if(read%10000 == 0)
            Message.info(read+ " lines read");
          Variant v = new Variant(line);
          if(v.isValid())
            variants.add(v);
        }
      }
      in.close();
      return variants;
    }
  }

  public static class Variant {
    private final String chr;
    private final int pos;
    private final String id;
    private final String[] alleles;
    private final String[] genotypes;

    public static final String MISSING = "-1";
    public static final String HETERO = "-2";

    public Variant(String line) {
      String[] f = line.split("\t");
      this.chr = f[0];
      this.pos = Integer.parseInt(f[1]);
      this.id = (f[2].isEmpty() || ".".equals(f[2])) ? chr+":"+pos : f[2];
      this.alleles = (f[3]+","+f[4]).split(",");
      this.genotypes = new String[f.length - 9];
      for(int i = 0; i < this.genotypes.length; i++)
        this.genotypes[i] = readGenotype(f[i+9]);
    }

    public boolean has(String chr, int pos, String allele, Mode mode){
      if(!this.chr.equalsIgnoreCase(chr))
        return false;
      if(this.pos != pos)
        return false;
      if(mode == Mode.IGNORE)
        return true;
      for(String al : alleles) {
        if (al.equalsIgnoreCase(allele) && mode == Mode.HOMOZYGOUS)
          return true;
        if (al.equalsIgnoreCase(HETERO) && mode == Mode.HETEROZYGOUS)
          return true;
      }
      return false;
    }

    private String readGenotype(String s) {
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

    public boolean canBeTarget(String targetAllele, Mode mode){
      for(String geno : genotypes) {
        if (!targetAllele.equalsIgnoreCase(geno) && mode == Mode.HOMOZYGOUS)
          return false;//all allele should be the same
        if (!targetAllele.equalsIgnoreCase(HETERO) && mode == Mode.HETEROZYGOUS)
          return false;//all allele should be the hetero
      }
      return true;
    }

    public boolean isValid(){
      for(String geno : genotypes){
        if(HETERO.equals(geno))
          return false;//allele should be homozygous
        if(MISSING.equals(geno))
          return false;//no missing allele
      }
      return true;
    }

    public String[] getGenotypes(){
      return genotypes;
    }

    public String getChr() {
      return chr;
    }

    public int getPos() {
      return pos;
    }

    public String getId() {
      return id;
    }

    public String[] getAlleles() {
      return alleles;
    }

    public void removeNonAncestral(ArrayList<Integer> samples){
      String noanc = "NO-ANCESTRAL";
      String anc = noanc;
      HashMap<String, Integer> counts = new HashMap<>();
      for(int i : samples){
        String allele = genotypes[i];
        Integer count = counts.get(allele);
        if(count == null)
          count = 0;
        counts.put(allele, ++count);
      }

      int max = -1;
      for(int count : counts.values()){
        if(count > max)
          max = count;
      }

      //Check if there is on max or many
      for(String allele : counts.keySet()){
        if(counts.get(allele) == max){
          if(anc.equals(noanc))//first
            anc = allele;
          else {//Multiple max, no ancestral
            samples.clear();
            return;
          }
        }
      }

      ArrayList<Integer> toRemove = new ArrayList<>();
      for(int i : samples){
        if(!anc.equals(genotypes[i]))
          toRemove.add(i);
      }

      samples.removeAll(toRemove);
    }
  }
}
