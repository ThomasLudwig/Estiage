package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * Class representing the human legible input file
 * 2 Versions:
 * RAW (just positions and genotypes)
 * COMPLETE (with recombination fractions and frequencies)
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2021-03-18
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class TSVFile {
  public enum Type{
    RAW,
    COMPLETE
  }

  /**
   * the marker we want to date
   */
  private Marker target;
  /**
   * the markers on the left (from closest to farthest to the marker to estimate)
   */
  private Marker[] leftMarkers;
  /**
   * the markers on the right (from closest to farthest to the marker to estimate)
   */
  private Marker[] rightMarkers;

  /**
   * List of Sample names
   */
  private final String[] samples;

  /**
   * Is format complete ? (Freq, theta mb)
   */
  private Type type;

  public TSVFile(String filename, Type type) throws IOException, EstiageFormatException {
    this.type = type;

    //read the whole file into an array
    ArrayList<String[]> tmp = new ArrayList<>();
    int nbCols = 0;
    UniversalReader in = new UniversalReader(filename);
    String line = in.readLine();
    String[] f = line.split("\t", -1);
    nbCols = f.length;
    tmp.add(f);
    int nbLines = 1;
    while ((line = in.readLine()) != null) {
      if(!line.isEmpty()) {
        nbLines++;
        f = line.split("\t", -1);
        if (f.length == nbCols)
          tmp.add(f);
        else
          throw new EstiageFormatException("Line [" + nbLines + "] has [" + f.length + "] columns, expected number [" + nbCols + "]");
      }
    }
    in.close();

    String[][] lines = new String[nbLines][nbCols];
    for(int i = 0 ; i < nbLines; i++) {
      String[] tt = tmp.get(i);
      for (int j = 0 ; j < nbCols; j++)
        lines[i][j] = tt[j];
    }

    //the whole file is in memory in a String[][]
    final int nbSamples = (type == Type.COMPLETE) ?
            nbLines - 10 :
            nbLines - 3;
    Message.info("Samples : "+nbSamples);
    final int lName = 1;
    final int lPos = 1 + nbSamples + 1;
    final int lAncestral = lPos + 1;
    final int lMb = lAncestral + 1;
    final int lDistance = lMb + 1;
    final int lMeanRate = lDistance + 1;
    final int lcM = lMeanRate + 1;
    final int lTheta = lcM + 1;
    final int lFreq = lTheta + 1;

    int nbLeft = 0;
    int nbRight = 0;

    try {
      nbLeft = Integer.parseInt(lines[0][1].toLowerCase().replace("left", ""));
    } catch(ArrayIndexOutOfBoundsException | NumberFormatException e){
      Message.fatal("Cell [0][1] should contain \"LeftX\" where X is the number of markers to the left of the target");
    }
    try{
      nbRight = Integer.parseInt(lines[0][nbCols - 1].toLowerCase().replace("right", ""));
    } catch(ArrayIndexOutOfBoundsException | NumberFormatException e){
      Message.fatal("Cell [0][COL-1] should contain \"RightX\" where X is the number of markers to the right of the target");
    }

    final int cTarget = nbLeft + 1;

    target = null;
    leftMarkers = new Marker[nbLeft];
    rightMarkers = new Marker[nbRight];
    samples = new String[nbSamples];

    for(int s = 0; s < nbSamples; s++)
      samples[s] = lines[s + 2][0];

    target = new Marker(lines[1][cTarget]);
    target.setChromosomeAndPosition(lines[lPos][cTarget]);
    if(type == Type.COMPLETE)
      target.setMegaBases(Double.parseDouble(lines[lMb][cTarget]));

    for(int c = 1 ; c < nbCols; c++){
      if(c == cTarget)
        continue;
      Marker m = new Marker(lines[lName][c]);
      m.setChromosomeAndPosition(lines[lPos][c]);
      if(type == Type.COMPLETE){
        m.setMegaBases(Double.parseDouble(lines[lMb][c]));
        m.setRecombinationFraction(Double.parseDouble(lines[lTheta][c]));
        m.setFrequencies(Double.parseDouble(lines[lFreq][c]));
      }
      for(int s = 0; s < nbSamples; s++)
        m.addHaplotype(samples[s], lines[2 + s][c]);

      m.setDistance(target);

      if(c < cTarget)
        leftMarkers[nbLeft - c] = m;
      else
        rightMarkers[c - (cTarget + 1)] = m;
    }

    //here set Ancestral Allele
    Message.info("Compute left ancestral");
    this.leftMarkers = this.computeAncestralAlleles(this.leftMarkers);
    Message.info("Compute right ancestral");
    this.rightMarkers = this.computeAncestralAlleles(this.rightMarkers);
  }

  public TSVFile(Marker target, Marker[] leftMarkers, Marker[] rightMarkers, String[] samples){
    this.target = target;
    this.leftMarkers = leftMarkers;
    this.rightMarkers = rightMarkers;
    this.samples = samples;
  }

  public static final String T = "\t";

  private Marker[] computeAncestralAlleles(Marker[] markers){
    //example
    //1 5 <-
    //1 5 <-
    //1 5 <-
    //1 5 <-
    //2 5 <-
    //2 5 <-
    //2 1 <-
    //2 1 <-
    //2 1 <-
    //2 1 <-
    //first keep only 5
    //1 5 <-
    //1 5 <-
    //1 5 <-
    //1 5 <-
    //2 5 <-
    //2 5 <-
    //the ancestry is 1 5 <-, not 2 5

    //compute ancestral allele A for marker m
    //For each sample that differs from A, set all further markers to missing
    for(int i = 0; i < markers.length; i++){
      Marker m = markers[i];
      m.setAncestral();
      String ancestral = m.getAncestral();
      String alleles = "";
      for(String val : m.getAllAlleles())
        alleles += ","+val;
      if(alleles.length() > 0)
        alleles = alleles.substring(1);
      Message.info("["+i+"] : ["+ancestral+"] {"+ alleles+"}");
      if("-1".equals(ancestral)) {
        Message.info("Last Marker ("+i+")");
        if(i == markers.length - 1)
          return markers;
        //Here we copy following markers if haplotypes have begun to diverge
        Marker[] ret = new Marker[i+1];
        System.arraycopy(markers, 0, ret, 0, i + 1);
        return ret;
      } else {
        ArrayList<String> diff = new ArrayList<>();
        String skip = "Skip :";
        for (String sample : samples) {
          String haplo = m.getHaplotype(sample);
          if (!haplo.isEmpty() && !haplo.equals(ancestral)) {
            diff.add(sample);
            skip += " " + sample;
          }
        }
        Message.info(skip);
        for (int j = i + 1; j < markers.length; j++)
          for (String sample : diff)
            markers[j].setHaplotype(sample, Marker.MISSING);
      }
    }
    return markers;
  }

  /**
   * Assigns Frequencies to all markers
   * @param gnomadFilename the VCF file containing the frequencies (in the INFO:AF annotation)
   * @param hapmapFilename the TSV file containing the mutation rate (in the 3rd column)
   * @throws IOException if there is a problem reading the files
   * @throws EstiageFormatException if there are missing data in the mutation rate Map
   */
  public void generateFrequenciesAndRecombinationFractions(String gnomadFilename, String hapmapFilename) throws IOException, EstiageFormatException {
    Marker first = target;
    Marker last = target;
    if(leftMarkers.length > 0)
      first = leftMarkers[leftMarkers.length - 1];
    if(rightMarkers.length > 0)
      last = rightMarkers[rightMarkers.length - 1];
    if(first != null && last != null){
      GnomAD gnomad = new GnomAD(gnomadFilename);
      Message.info("Reading mutation rates");
      HapMap hapmap = new HapMap(hapmapFilename, first.getPosition(), last.getPosition());
      Message.info("Applying");
      for(Marker m : leftMarkers) {
        gnomad.applyFrequency(m);
        hapmap.applyRate(m, target);
      }
      for(Marker m : rightMarkers) {
        gnomad.applyFrequency(m);
        hapmap.applyRate(m, target);
      }
    }
    type = Type.COMPLETE;
  }

  public Marker[] getLeftMarkers() { return leftMarkers; }

  public Marker[] getRightMarkers() {
    return rightMarkers;
  }

  public String[] getSamples() {
    return samples;
  }

  public void export(String filename) throws IOException {
    PrintWriter out = new PrintWriter(new FileWriter(filename));
    String header  = ".";
    String nameLine = "Samples\\Markers";
    String[] samplesLines = (String[])samples.clone();
    String positionLine="Position";
    String ancestralLine="Ancestral";
    String mbLine="mb";
    String distanceLine="Distance";
    String rateLine="Mean rate";
    String cMLine="cM";
    String thetaLine="θ Recombination Fraction";
    String freqLine="Freq";
    Marker m;

    for(int i = leftMarkers.length-1; i >= 0; i--) {
      header += T + "Left" + (i + 1);
      m = leftMarkers[i];
      nameLine += T + m.getName();

      for(int s = 0; s < samples.length; s++)
        samplesLines[s] += T + m.getHaplotype(samples[s]);
      positionLine += T + m.getChromosome()+":"+m.getPosition();
      if (type == Type.COMPLETE) {
        ancestralLine += T + m.getAncestral();
        mbLine += T + m.getMegaBases();
        distanceLine += T + m.getDistance();
        rateLine += T + m.getRate();
        cMLine += T + m.getcM();
        thetaLine += T + m.getRecombinationFraction();
        freqLine += T + m.getFrequencies();
      }
    }

    m = target;
    header+= T + "Target";
    nameLine+= T + m.getName();
    for(int s = 0; s < samples.length; s++)
      samplesLines[s] += T;
    positionLine += T + m.getChromosome()+":"+m.getPosition();
    if(type == Type.COMPLETE) {
      ancestralLine += T;
      mbLine += T + m.getMegaBases();
      distanceLine += T;
      rateLine += T;
      cMLine += T;
      thetaLine += T;
      freqLine += T;
    }

    for(int i = 0; i < rightMarkers.length; i++){
      header += T + "Right"+(i+1);
      m = rightMarkers[i];
      nameLine += T + m.getName();
      for(int s = 0; s < samples.length; s++)
        samplesLines[s] += T + m.getHaplotype(samples[s]);
      positionLine += T + m.getChromosome()+":"+m.getPosition();
      if(type == Type.COMPLETE){
        ancestralLine += T + m.getAncestral();
        mbLine += T + m.getMegaBases();
        distanceLine += T + m.getDistance();
        rateLine += T + m.getRate();
        cMLine += T + m.getcM();
        thetaLine += T + m.getRecombinationFraction();
        freqLine += T + m.getFrequencies();
      }
    }

    out.println(header);
    out.println(nameLine);
    for(String sampleLine : samplesLines)
      out.println(sampleLine);
    out.println(positionLine);

    if(type == Type.COMPLETE){
      out.println(ancestralLine);
      out.println(mbLine);
      out.println(distanceLine);
      out.println(rateLine);
      out.println(cMLine);
      out.println(thetaLine);
      out.println(freqLine);
    }

    out.close();
  }

  public void printSummary(){
    String message = "Target "+target.getName()+" ("+target.getChromosome()+":"+target.getPosition()+")";
    Marker left = leftMarkers[leftMarkers.length-1];
    Marker right = rightMarkers[rightMarkers.length-1];
    message += "\nLeft("+leftMarkers.length+") "+left.getName()+" ("+left.getChromosome()+":"+left.getPosition()+")";
    message += "\nRight("+rightMarkers.length+") "+right.getName()+" ("+right.getChromosome()+":"+right.getPosition()+")";
    Message.info(message);
  }
}