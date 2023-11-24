package fr.inserm.u1078.estiage;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Ranks elements (sorted be number of occurrences DESC)
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-23
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Ranking<T> {
  //TODO generalize sorting order and put in a lib

  private final Object[] rankedElements;
  private final int[] occurrences;

  /**
   * Builds a Ranking from a HashMap(Element, Number of occurrences)
   * @param counts keys are elements, values are number of occurrences for each element
   */
  public Ranking(HashMap<T, Integer> counts) {
    int size = counts.size();
    rankedElements = new Object[size];
    occurrences = new int[size];

    //rank the values by count DESC
    for (int i = 0; i < size; i++) {
      int maxOccurrence = 0;
      T top = null;
      for (T element : counts.keySet()) {
        int occurrence = counts.get(element);
        if (occurrence > maxOccurrence) {
          maxOccurrence = occurrence;
          top = element;
        }
      }
      occurrences[i] = maxOccurrence;
      rankedElements[i] = top;
      counts.remove(top);
    }
  }

  /**
   * Gets the value of the ith element
   * @param i the index (rank) of the element
   * @return the element
   */
  @SuppressWarnings("unchecked")
  public T getElement(int i) {
    return (T) rankedElements[i];
  }

  /**
   * Gets the number of occurrences for the ith element
   * @param i the index (rank) of the element
   * @return the number of occurrences
   */
  public int getOccurrence(int i) {
    return occurrences[i];
  }

  /**
   * Gets the size of the ranking
   * @return the number of ranked elements
   */
  public int size(){
    return rankedElements.length;
  }

  /**
   * Gets all elements with top occurrences
   * @return an ArrayList with all the top elements
   */
  public ArrayList<T> getTopElements(){
    ArrayList<T> ret = new ArrayList<>();
    if(!isEmpty()) {
      ret.add(getElement(0));
      for(int i = 1 ; i < size(); i++)
        if(getOccurrence(i) == getOccurrence(0))
          ret.add(getElement(i));
        else
          break;
    }
    return ret;
  }

  /**
   * Checks if there are multiple top elements
   * @return true if there is more than 1 element with max occurrences
   */
  public boolean hasTopExAequo(){
    if(size() < 2)
      return false;
    return occurrences[0] == occurrences[1];
  }

  /**
   * Check if the ranking is empty
   * @return true if the ranking is empty
   */
  public boolean isEmpty(){
    return size() < 1;
  }
}
