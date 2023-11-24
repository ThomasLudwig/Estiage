package fr.inserm.u1078.estiage.ctranslation;

/**
 * Exception during the call to Estiage
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-20
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class EstiageException extends Exception {
  public EstiageException(String message) {
    super(message);
  }
}
