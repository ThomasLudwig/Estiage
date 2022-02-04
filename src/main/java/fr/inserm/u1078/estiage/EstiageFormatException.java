package fr.inserm.u1078.estiage;

/**
 * Exception Thrown when data are not formatted as expected
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2021-03-18
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class EstiageFormatException extends Exception {

  public EstiageFormatException() {
  }

  public EstiageFormatException(String message) {
    super(message);
  }

  public EstiageFormatException(String message, Throwable cause) {
    super(message, cause);
  }

  public EstiageFormatException(Throwable cause) {
    super(cause);
  }

  public EstiageFormatException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
    super(message, cause, enableSuppression, writableStackTrace);
  }
}
