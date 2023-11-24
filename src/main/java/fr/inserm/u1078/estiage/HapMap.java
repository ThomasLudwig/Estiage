package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.IOException;
import java.util.*;

/**
 * Class used to request mutation rates from HapMap
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2022-01-11
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class HapMap {
  private final TreeMap<Integer, Double> mutationRates;

  /**
   * Creates an HapMap object and loads the data
   * @param filename the name of the HapMap file
   * @param first the first position to consider
   * @param last the last position to consider
   * @throws IOException
   */
  public HapMap(String filename, int first, int last) throws IOException {
    this.mutationRates = HapMap.load(filename, first, last);
    if(mutationRates == null){
      Message.error("Mutation Rates is null after load");
    }
  }

  /**
   * Load the HapMap mutation rates
   * The chromosome is not needed, as the files are expected to be split by chr
   * @param filename the name of the the hapmap file
   * @param first the position of the first marker
   * @param last the position of the last marker
   * @return the loaded TreeMap of mutation rates
   * @throws IOException if the file can't be read
   */
  private static TreeMap<Integer, Double> load(String filename, int first, int last) throws IOException {
    Message.info("Looking in ["+filename+"] from "+filename+" between "+first+"-"+last);
    UniversalReader in = new UniversalReader(filename);
    TreeMap<Integer, Double> mutationRates = new TreeMap<>();
    String line = in.readLine(); // skip header
    /*
    Format
    Chromosome  Position(bp)  Rate(cM/Mb) Map(cM)
    chr17       13043         3.745829    0.000000
    chr17       41084         3.748065    0.105037
    chr17       51088         2.578174    0.142532
    chr17       52467         2.573723    0.146088
    chr17       53011         2.532886    0.147488
    chr17       53206         2.487219    0.147982
    chr17       57026         2.494051    0.157483
    chr17       58275         2.496189    0.160598
    chr17       58449         2.500471    0.161032

     */

    //Need to read line BEFORE first marker (hard, buffer each line and unbuffer after first marker) and Line AFTER last marker (easy)
    //Gradiant rate between each position (P1-P2)/N
    //Mean from gradiant between each marker
    int prevPos = -1;
    double prevRate = -1D;
    int read = 0;
    long start = new Date().getTime();
    boolean hasF = false;
    while((line = in.readLine()) != null) {
      read++;
      if(read%100000 == 0){
        double dur = (new Date().getTime() - start)/1000D;
        int speed = (int)(read/dur);
        Message.info("Read : "+read+" lines in "+((int)dur)+" s. "+speed+" l/s");
      }
      String[] f = line.split("\t", -1);
      try {
        int pos = Integer.parseInt(f[1]);
        double rate = Double.parseDouble(f[2]);
        if (pos > first) { //everytime after first
          if(!hasF) {
            Message.info("Start because " + pos + " > " + first);
            hasF = true;
          }
          mutationRates.put(prevPos, prevRate); //put previous
        }
        if (pos > last) { //on last value, override buffer
          mutationRates.put(pos, rate); //add last and quit
          Message.info("Stop because "+pos+" > "+last);
          break;
        }
        //update previous
        prevPos = pos;
        prevRate = rate;
      } catch(Exception e){
        Message.error("could not parse line ["+line+"]", e);
      }
    }

    //Loading done
    if(mutationRates.isEmpty())
      Message.error("Empty rate list");
    int lowest = mutationRates.firstKey();
    int highest = mutationRates.lastKey();
    Message.info("First : ["+lowest+" / "+mutationRates.get(lowest)+"]");
    Message.info("Last : ["+highest+" / "+mutationRates.get(highest)+"]");

    return mutationRates;
  }

  /**
   * Add a position for an existing marker. If the position already exist, nothing to do. Otherwise, compute the rate at the position
   * @param pos
   * @throws EstiageFormatException
   */
  private void addMarker(int pos) throws EstiageFormatException {
    if(mutationRates.isEmpty())
      throw new EstiageFormatException("Hapmap data are empty");
    if(mutationRates.keySet().contains(pos))//Value added
      return;

    int first = mutationRates.firstKey();
    int last = mutationRates.lastKey();

    if(pos < first) {
      double v = mutationRates.get(first);
      mutationRates.put(pos, v);
      return ;
    }

    if(last < pos) {
      double v = mutationRates.get(last);
      mutationRates.put(pos, v);
      return;
    }

    Iterator<Integer> keys = mutationRates.navigableKeySet().iterator();

    int left = keys.next();
    while(keys.hasNext()){
      int right = keys.next();
      if(left < pos && pos < right) {
        double leftValue = mutationRates.get(left);
        double rightValue = mutationRates.get(right);
        int interval = right - left; //can't be 0 since left < pos && pos < right, so left != right
        double gradiant = (rightValue - leftValue) / interval;
        double value = leftValue + gradiant * (pos - left);
        mutationRates.put(pos, value);
        return;
      }
      left = right;
    }

    throw new EstiageFormatException("Could not add position ["+pos+"] in hapmap data ["+first+";"+last+"]");
  }

  /**
   * Gets the list of available positions in [p1;p2]
   * @param p1 the first position
   * @param p2 the second position
   * @return
   */
  private ArrayList<Integer> getPositions(int p1, int p2){
    ArrayList<Integer> ret = new ArrayList<>();
    for(Integer pos : mutationRates.navigableKeySet())
      if(p1 <= pos && pos <= p2)
        ret.add(pos);
    return ret;
  }

  /**
   * Gets the mean recombination rate between two unsorted positions
   * @param p1 the first position
   * @param p2 the second position
   * @return the recombination rate
   * @throws EstiageFormatException
   */
  public double getRate(int p1, int p2) throws EstiageFormatException {
    //Message.debug("From ["+p1+"] to ["+p2+"]");
    //add points if missing
    addMarker(p1);
    addMarker(p2);
    //Message.debug("At["+p1+"]>["+mutationRates.get(p1)+"], At["+p2+"]>["+mutationRates.get(p2)+"]");
    //get all points in interval
    int first = Math.min(p1, p2);
    int last = Math.max(p1, p2);
    ArrayList<Integer> points = getPositions(first, last);

    double rateSum = 0;

    double rightValue = 0;
    for(int i = 1; i < points.size(); i++){
      int left = points.get(i-1);
      int right = points.get(i);
      int d = right - left; //d != 0 because right > left, as points are unique
      double leftValue = mutationRates.get(left);
      rightValue = mutationRates.get(right);

      double mean = 0.5 * (rightValue + leftValue);
      //Message.debug("["+left+" ; "+right+"[ >  ["+mean+"]");
      rateSum += d * mean; //here, add [left;right[, to avoid adding twice each middle point
    }
    //here, add right
    rateSum += rightValue;
    int distance = 1 + last - first;
    //Message.debug("Mean["+(rateSum / distance)+"]=["+rateSum+"]/["+distance+"]");
    return rateSum / distance;

  }

  /**
   * Gets the mean recombination rate between two unsorted markers
   * @param m1 the first marker
   * @param m2 the second marker
   * @return the mean recombination rate
   * @throws EstiageFormatException
   */
  public double getRate(Marker m1, Marker m2) throws EstiageFormatException {
    return getRate(m1.getPosition(), m2.getPosition());
  }
}
