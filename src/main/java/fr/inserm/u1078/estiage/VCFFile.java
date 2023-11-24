package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.*;
import java.util.ArrayList;
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
  public enum Method {CLASSICAL, LONGEST_HAPLOTYPE}

  private final String filename;
  private final boolean isTabix;
  private final String[] samples;
  private final Mode mode;

  private String chr;
  private int pos;
  private String allele;

  /**
   * Creates a VCFFile Object from a VCF file
   * @param filename the name of the VCF file
   * @param mode the mode to use
   * @throws IOException the VCF file can't be read
   * @throws EstiageFormatException the VCF file has no header or not enough columns
   */
  public VCFFile(String filename, Mode mode) throws IOException, EstiageFormatException {
    this.filename = filename;
    this.mode = mode;
    this.isTabix = isTabix();
    this.samples = readSamples();
  }

  /**
   * Does the tabix file exist ?
   * @return true, if the tabix file exists
   * @throws FileNotFoundException if the VCF file doesn't exist
   */
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

  /**
   * Gets the list of samples from the header
   * @return all the samples in an array of String
   * @throws IOException the VCF file can't be read
   * @throws EstiageFormatException the VCF file has no header or not enough columns
   */
  private String[] readSamples() throws IOException, EstiageFormatException {
    UniversalReader in = new UniversalReader(this.filename);
    String line;
    String header = null;
    while((line = in.readLine()) != null){
      if(line.startsWith("#")) {
        if(line.startsWith("#CHROM"))
          header = line;
        //else, skip irrelevant headers
      }
      else
        break;
    }
    in.close();
    if(header == null)
      throw new EstiageFormatException("No header found in VCF ["+filename+"]");
    String[] f = header.split("\t");
    if(f.length < 9)
      throw new EstiageFormatException("The VCF file ["+filename+"] only has ["+f.length+"] columns");
    String[] s = new String[f.length - 9];
    System.arraycopy(f, 9, s, 0, s.length);
    return s;
  }

  /**
   * Set the chromosome, position and allele of the variant of interest
   * @param chrPosAllele - a String in the format "chrom:position:allele"
   */
  public void setVariant(String chrPosAllele) {
    String[] f = chrPosAllele.split(":");
    this.chr = f[0];
    this.pos = Integer.parseInt(f[1]);
    this.allele = f[2];
  }

  /**
   * Gets the indices of the first left/right variants. Here "first" means closest to the variant of interest. If starts looking from the middle of the array
   * @param variants the valid variants for the chromosome
   * @return an integer array with the values {indexLeft, indexRight}
   */
  private int[] getLeftRight(ArrayList<VCFVariant> variants){
    int size = variants.size();
    int left = Integer.MAX_VALUE;
    int right = -1;

    int start = 0;
    int end = variants.size() - 1;

    while(start <= end ) {
      int i = (start + end) / 2;
      //Message.info("["+start+"|"+i+"|"+end+"]");
      int v = variants.get(i).getPos();
      if (v == pos) {
        return new int[]{i - 1, i + 1};
      } else if (v < pos) {
        if (i == size - 1)
          return new int[]{i, i + 1};
        if (variants.get(i + 1).getPos() > pos)
          return new int[]{i, i + 1};
        start = i;
      } else { //v > pos
        if (i == 0)
          return new int[]{i - 1, i};
        if (variants.get(i - 1).getPos() < pos)
          return new int[]{i - 1, i};
        end = i;
      }
    }

    return new int[]{left, right};
  }

  /**
   * The "main" method. It reads the VCF file and exports it as a raw file
   * @param raw the name of the output file
   * @param method CLASSICAL or LONGEST_HAPLOTYPE
   * @throws IOException if the VCF file can't be read
   * @throws EstiageFormatException the VCF file has no header or not enough columns, if the VCF File is unsorted, if the genotype for a variant doesn't match the selected mode
   * @throws InterruptedException if there is a problem with the thread calling tabix
   */
  public void exportAsRaw(String raw, Method method) throws IOException, EstiageFormatException, InterruptedException {
    ArrayList<VCFVariant> variants = loadChromosome();
    Message.info(variants.size() + " valid variants found in [" + this.filename + "] on chromosome [" + this.chr + "]");
    int[] leftRight = getLeftRight(variants);
    int left = leftRight[0];
    int right = leftRight[1];
    //Message.info("Looking for ["+pos+"] and found here ["+left+","+right+"] ["+variants.get(left).pos+";"+variants.get(right).pos+"]");
    if (left > right)
      throw new EstiageFormatException("Something went wrong. Is you VCF file unsorted?");
    Marker target;
    if (left + 2 == right) {
      VCFVariant targetVariant = variants.get(left + 1);
      //Message.info("Target is "+targetVariant.pos);
      if (!targetVariant.has(chr, pos, allele, mode))
        throw new EstiageFormatException("The genotypes for the variant at [" + chr + ":" + pos + ":" + allele + "] do not match the selected mode [" + mode + "]");
      target = new Marker(targetVariant, samples);
    } else {
      target = new Marker(chr + ":" + pos + ":" + allele, samples);
      //Message.info("Target is default");
    }

    //Here we have the Target Variant, and the index of the first lefT/right variants
    Message.info("Target Variant is between  [" + (left + 1) + ";" + (right + 1) + "]/" + variants.size());
    //List of samples not yet excluded
    ArrayList<Integer> leftSamples = new ArrayList<>();
    ArrayList<Integer> rightSamples = new ArrayList<>();
    for (int s = 0; s < samples.length; s++) {
      leftSamples.add(s);
      rightSamples.add(s);
    }

    //Left
    //Rebuilding the longest left side possible for the haplotype
    TreeMap<Integer, VCFVariant> leftVariants = new TreeMap<>();
    for(int i = left; i >= 0 && !leftSamples.isEmpty(); i--)
      processVariant(variants.get(i), leftSamples, leftVariants, method);
    //Builds Markers from LeftVariant + Samples
    ArrayList<Marker> leftMarkers = buildMarkers(leftVariants);
    Message.info("Variants kept on the left side : "+leftMarkers.size());

    //Right
    //Rebuilding the longest right side possible for the haplotype
    TreeMap<Integer, VCFVariant> rightVariants = new TreeMap<>();
    for(int i = right ; i < variants.size() && !rightSamples.isEmpty(); i++)
      processVariant(variants.get(i), rightSamples, rightVariants, method);
    //Builds Markers from RightVariant + Samples
    ArrayList<Marker> rightMarkers = buildMarkers(rightVariants);
    Message.info("Variants kept on the right side : "+rightMarkers.size());

    //Export TSV
    TSVFile tsvFile = new TSVFile(target, leftMarkers.toArray(new Marker[0]), rightMarkers.toArray(new Marker[0]), samples);
    tsvFile.export(raw);
  }

  /**
   * For the current variants, removes the samples that do not have the ancestralAllele and update the variant list
   * @param v the current variant
   * @param remainingSamples the list of remainingSamples (will be updated)
   * @param sideVariants the TreeMap of Variants, the key is the number of remaining samples
   * @param method CLASSIC or LONGEST_HAPLOTYPE
   */
  private void processVariant(VCFVariant v, ArrayList<Integer> remainingSamples, TreeMap<Integer, VCFVariant> sideVariants, Method method) {
    //Get the topAlleles for the Variant
    String[] topAlleles = v.getTopAllelesAndCount(remainingSamples);
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++
      + Start in the difference with the classical method +
      +++++++++++++++++++++++++++++++++++++++++++++++++++++*/
    /*
    When following the haplotype, for variants V
    1. There is a majority (ancestral variant)            -> stop haplotype for other samples (add -1 of letter)
    2. There is a tie, and each top allele has count == 1 -> stop haplotype for all samples   (add -1 of letter)
    3. There is a tie, and top allele has count > 1       -> ignore this variant
     */
    boolean drop = topAlleles.length > 2 && Integer.parseInt(topAlleles[0]) > 1;
    if(drop && method == Method.LONGEST_HAPLOTYPE){
      Message.debug("Drop the variant ["+v.getChr()+":"+v.getPos()+"] number of topAllele ("+(topAlleles.length-1)+") for {"+topAlleles[0]+"} samples");
    } else {
      //Update the list of remaining samples
      v.removeNonAncestral(remainingSamples, topAlleles);
      //put the current variant in the list of markers
      int size = remainingSamples.size();
      sideVariants.put(size, v);
    }
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++
      + End in the difference with the classical method +
      +++++++++++++++++++++++++++++++++++++++++++++++++++*/
  }

  //WHY keep FIRST and LAST, why not just keep FIRST ? or LAST ?
  //If think the LAST would be enough and to would make a simpler code (no TreeMap) and smaller TSV (less columns)
  // After some tests, it works better with only the LAST variant

  /**
   * Migrates the Variants from the TreeMap to the Marker List.<br/>
   * The intermediate TreeMap is to take the FIRST and LAST variant of each "size" and drop those in the middle
   * @param variants a TreeMap of Variant List, where the number of samples with the ancestral allele is the key
   * @return the ordered list of Markers from the Variants
   */
  private ArrayList<Marker> buildMarkers(TreeMap<Integer, VCFVariant> variants) {
    ArrayList<Marker> markers = new ArrayList<>();
    for(int key : variants.descendingKeySet())
      markers.add(new Marker(variants.get(key), samples));

    return markers;
  }

  /**
   * Gets a list of variants from a list of VCF lines
   * @param vcfLines a list of VCF lines
   * @param onlyValid if true, only valid variants are kept
   * @return the list of variants
   * @throws IOException if the VCF file can't be read
   * @throws InterruptedException if there is a problem with the thread calling tabix
   */
  private ArrayList<VCFVariant> getVariants(ArrayList<String> vcfLines, boolean onlyValid) throws IOException, InterruptedException {
    ArrayList<VCFVariant> variants = new ArrayList<>();
    for(String vcfLine : vcfLines){
      VCFVariant v = new VCFVariant(vcfLine);
      if(!onlyValid || v.isValid())
        variants.add(v);
    }
    return variants;
  }

  /**
   * If all is well, does nothing. Throws exceptions or prints warning otherwise
   * @throws EstiageFormatException if the variant is not in the file of if the genotype doesn't match the selected mode
   * @throws IOException if the VCF file can't be read
   * @throws InterruptedException if there is a problem with the thread calling tabix
   */
  private void checkHasVariant() throws EstiageFormatException, IOException, InterruptedException {
    if(this.mode == Mode.IGNORE)
      return;
    if(!this.isTabix) {
      Message.warning("As the file is not tabixed, no check will be performed on the target variant's genotypes");
      return;
    }

    for(VCFVariant v : this.getVariants(Utils.getLinesFromTabixedVCF(filename, chr, pos), false)) {
      if (v.has(chr, pos, allele, mode)) {
        if (!v.canBeTarget(allele, mode))
          throw new EstiageFormatException("Variant to estimate should be homozygous for each sample in the file. Here: " + String.join(",", v.getGenotypes()));
        return; //Found
      }
    }
    throw new EstiageFormatException("VCF file [" + filename + "] doesn't not contains a variant for " + chr + ":" + pos + ":" + allele);
  }

  /**
   * Loads of the valid variants for the chromosome
   * @return ArrayList of variants
   * @throws EstiageFormatException if the variant is not in the file of if the genotype doesn't match the selected mode
   * @throws IOException if the VCF file can't be read
   * @throws InterruptedException if there is a problem with the thread calling tabix
   */
  private ArrayList<VCFVariant> loadChromosome() throws IOException, EstiageFormatException, InterruptedException {
    this.checkHasVariant();
    if(this.isTabix){
      return this.getVariants(Utils.getLinesFromTabixedVCF(filename, chr), true);
    } else {
      ArrayList<VCFVariant> variants = new ArrayList<>();
      Message.warning("File ["+filename+"] is not tabixed, this will be slow");
      UniversalReader in = new UniversalReader(this.filename);
      String line;
      int read = 0;
      while((line = in.readLine()) != null){
        if(line.startsWith(chr)){
          read++;
          if(read%10000 == 0)
            Message.info(read+ " lines read");
          VCFVariant v = new VCFVariant(line);
          if(v.isValid())
            variants.add(v);
        }
      }
      in.close();
      return variants;
    }
  }
}
