package fr.inserm.u1078.estiage.test;

import fr.inserm.u1078.estiage.Main;
import fr.inserm.u1078.tludwig.maok.UniversalReader;
import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;

/**
 * XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2023-11-21
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class TestAlgorithm {
  public static void main(String[] args) throws Exception {
    Message.setDebugActive(true);
    testPhasing();
    //testF508Del();
  }

  private static void testPhasing() throws Exception {
    String inputData = "Marker\tD1S211\t\tD1S451\t\tD1S2720\t\tD1S197\t\tD1S2661\t\tD1S417\t\trs28942111\t\tD1S200\t\tD1S2742\t\tD1S220\t\tD1S473\t\tD1S390\t\n" +
            "HC2\t163\t185\t175\t177\t237\t241\t132\t136\t182\t186\t188\t190\tA\tA\t169\t171\t242\t248\t235\t237\t236\t238\t196\t204\n" +
            "HC806\t169\t183\t175\t175\t241\t241\t140\t140\t\t\t188\t188\tA\tA\t169\t169\t248\t248\t235\t235\t240\t242\t204\t208\n" +
            "HC92\t181\t183\t175\t177\t235\t235\t134\t134\t182\t182\t188\t190\tA\tA\t169\t169\t248\t248\t223\t223\t242\t242\t200\t204\n" +
            "HC748\t183\t183\t175\t175\t237\t237\t138\t138\t182\t182\t188\t188\tA\tA\t169\t169\t260\t260\t235\t235\t240\t240\t208\t208\n" +
            "HC2062\t165\t165\t175\t177\t235\t237\t132\t136\t\t\t188\t188\tA\tA\t169\t169\t248\t248\t235\t235\t240\t240\t208\t208\n" +
            "Lyon1\t175\t177\t175\t177\t237\t241\t132\t132\t186\t186\t188\t190\tA\tA\t169\t169\t246\t246\t221\t225\t252\t252\t200\t204\n" +
            "Lyon2\t163\t163\t175\t175\t237\t237\t142\t142\t186\t188\t188\t188\tA\tA\t169\t169\t246\t246\t223\t233\t236\t236\t212\t212\n" +
            "S Afr\t165\t175\t175\t175\t237\t237\t138\t138\t182\t182\t188\t188\tA\tA\t167\t167\t256\t258\t237\t237\t242\t242\t200\t200\n" +
            "CAD3077\t163\t163\t175\t175\t237\t237\t136\t136\t186\t186\t188\t188\tA\tT\t175\t175\t248\t248\t233\t233\t236\t236\t212\t212\n" +
            "CAD3087\t165\t165\t177\t177\t237\t237\t132\t132\t190\t190\t190\t190\tA\tT\t169\t169\t246\t246\t235\t235\t236\t236\t208\t208\n" +
            "CAD3242\t183\t183\t177\t177\t237\t237\t136\t136\t190\t190\t190\t190\tA\tT\t173\t173\t244\t244\t223\t223\t242\t242\t204\t204\n" +
            "CAD3428\t163\t163\t175\t175\t237\t237\t134\t134\t186\t186\t188\t188\tA\tT\t169\t169\t246\t246\t223\t223\t242\t242\t208\t208\n" +
            "CAD3553\t183\t183\t179\t179\t237\t237\t134\t134\t186\t186\t188\t188\tA\tT\t169\t169\t248\t248\t237\t237\t244\t244\t204\t204\n" +
            "CAD3645\t163\t163\t175\t175\t237\t237\t138\t138\t182\t182\t188\t188\tA\tT\t169\t169\t248\t248\t225\t225\t242\t242\t208\t208\n" +
            "CAD3898\t183\t183\t175\t175\t237\t237\t134\t134\t188\t188\t188\t188\tA\tT\t169\t169\t248\t248\t223\t223\t244\t244\t196\t196\n" +
            "CAD4978\t185\t185\t179\t179\t237\t237\t142\t142\t188\t188\t188\t188\tA\tT\t169\t169\t248\t248\t235\t235\t244\t244\t208\t208\n" +
            "8053\t183\t183\t175\t175\t237\t237\t136\t136\t184\t184\t188\t188\tA\tT\t169\t169\t246\t246\t235\t235\t240\t240\t212\t212\n" +
            "3908\t165\t165\t173\t173\t237\t237\t128\t128\t182\t182\t188\t188\tA\tT\t163\t163\t240\t240\t235\t235\t244\t244\t196\t196\n" +
            "L2E05\t205\t205\t177\t177\t235\t235\t142\t142\t188\t188\t188\t188\tA\tT\t169\t169\t250\t250\t233\t233\t244\t244\t204\t204\n" +
            "CAD5914\t183\t183\t167\t167\t235\t235\t144\t144\t188\t188\t188\t188\tA\tT\t169\t169\t258\t258\t235\t235\t244\t244\t204\t204\n" +
            "NW1\t169\t169\t177\t177\t237\t237\t\t\t186\t186\t188\t188\tA\tT\t161\t161\t244\t244\t225\t225\t236\t236\t204\t204\n" +
            "NW3\t169\t169\t177\t177\t237\t237\t138\t138\t186\t186\t188\t188\tA\tT\t161\t161\t244\t244\t235\t235\t236\t236\t204\t204";
    final String expectedResultsTF = ".\tLeft3\tLeft2\tLeft1\tTarget\tRight1\tRight2\tRight3\tRight4\tRight5\n" +
            "Samples\\Markers\tD1S197\tD1S2661\tD1S417\trs28942111\tD1S200\tD1S2742\tD1S220\tD1S473\tD1S390\n" +
            "HC2\t132/136\t186\t188\tA\t169\t248\t235\t236/238\t\n" +
            "HC806\t\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "HC92\t\t182\t188\tA\t169\t248\t223\t\t\n" +
            "HC748\t\t182\t188\tA\t169\t260\t\t\t\n" +
            "HC2062\t\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "Lyon1\t132\t186\t188\tA\t169\t246\t\t\t\n" +
            "Lyon2\t142\t186\t188\tA\t169\t246\t\t\t\n" +
            "S Afr\t\t182\t188\tA\t167\t\t\t\t\n" +
            "CAD3077\t136\t186\t188\tA\t175\t\t\t\t\n" +
            "CAD3087\t\t\t190\tA\t169\t246\t\t\t\n" +
            "CAD3242\t\t\t190\tA\t173\t\t\t\t\n" +
            "CAD3428\t134\t186\t188\tA\t169\t246\t\t\t\n" +
            "CAD3553\t134\t186\t188\tA\t169\t248\t237\t\t\n" +
            "CAD3645\t\t182\t188\tA\t169\t248\t225\t\t\n" +
            "CAD3898\t\t188\t188\tA\t169\t248\t223\t\t\n" +
            "CAD4978\t\t188\t188\tA\t169\t248\t235\t244\t\n" +
            "8053\t\t184\t188\tA\t169\t246\t\t\t\n" +
            "3908\t\t182\t188\tA\t163\t\t\t\t\n" +
            "L2E05\t\t188\t188\tA\t169\t250\t\t\t\n" +
            "CAD5914\t\t188\t188\tA\t169\t258\t\t\t\n" +
            "NW1\t0\t186\t188\tA\t161\t\t\t\t\n" +
            "NW3\t138\t186\t188\tA\t161\t\t\t\t\n" +
            "Positions\t???:???\t???:???\t???:???\t55044016\t???:???\t???:???\t???:???\t???:???\t???:???";
    final String expectedResultsTT = ".\tLeft3\tLeft2\tLeft1\tTarget\tRight1\tRight2\tRight3\tRight4\tRight5\n" +
            "Samples\\Markers\tD1S197\tD1S2661\tD1S417\trs28942111\tD1S200\tD1S2742\tD1S220\tD1S473\tD1S390\n" +
            "HC2\t132/136\t186\t188\tA\t169\t248\t235\t236/238\t\n" +
            "HC806\t140\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "HC92\t\t182\t188\tA\t169\t248\t223\t\t\n" +
            "HC748\t\t182\t188\tA\t169\t260\t\t\t\n" +
            "HC2062\t132/136\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "Lyon1\t132\t186\t188\tA\t169\t246\t\t\t\n" +
            "Lyon2\t142\t186\t188\tA\t169\t246\t\t\t\n" +
            "S Afr\t\t182\t188\tA\t167\t\t\t\t\n" +
            "CAD3077\t136\t186\t188\tA\t175\t\t\t\t\n" +
            "CAD3087\t\t\t190\tA\t169\t246\t\t\t\n" +
            "CAD3242\t\t\t190\tA\t173\t\t\t\t\n" +
            "CAD3428\t134\t186\t188\tA\t169\t246\t\t\t\n" +
            "CAD3553\t134\t186\t188\tA\t169\t248\t237\t\t\n" +
            "CAD3645\t\t182\t188\tA\t169\t248\t225\t\t\n" +
            "CAD3898\t\t188\t188\tA\t169\t248\t223\t\t\n" +
            "CAD4978\t\t188\t188\tA\t169\t248\t235\t244\t\n" +
            "8053\t\t184\t188\tA\t169\t246\t\t\t\n" +
            "3908\t\t182\t188\tA\t163\t\t\t\t\n" +
            "L2E05\t\t188\t188\tA\t169\t250\t\t\t\n" +
            "CAD5914\t\t188\t188\tA\t169\t258\t\t\t\n" +
            "NW1\t0\t186\t188\tA\t161\t\t\t\t\n" +
            "NW3\t138\t186\t188\tA\t161\t\t\t\t\n" +
            "Positions\t???:???\t???:???\t???:???\t55044016\t???:???\t???:???\t???:???\t???:???\t???:???";
    String expectedResultsFT = ".\tLeft6\tLeft5\tLeft4\tLeft3\tLeft2\tLeft1\tTarget\tRight1\tRight2\tRight3\tRight4\tRight5\n" +
            "Samples\\Markers\tD1S211\tD1S451\tD1S2720\tD1S197\tD1S2661\tD1S417\trs28942111\tD1S200\tD1S2742\tD1S220\tD1S473\tD1S390\n" +
            "HC2\t163\t175/177\t237\t132/136\t186\t188\tA\t169\t248\t235\t236/238\t\n" +
            "HC806\t\t\t\t140\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "HC92\t\t\t\t\t182\t188\tA\t169\t248\t223\t\t\n" +
            "HC748\t\t\t\t\t182\t188\tA\t169\t260\t\t\t\n" +
            "HC2062\t165\t175/177\t237\t132/136\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "Lyon1\t175/177\t175/177\t237\t132\t186\t188\tA\t169\t246\t\t\t\n" +
            "Lyon2\t\t\t\t142\t186\t188\tA\t169\t246\t\t\t\n" +
            "S Afr\t\t\t\t\t182\t188\tA\t167\t\t\t\t\n" +
            "CAD3077\t163\t175\t237\t136\t186\t188\tA\t175\t\t\t\t\n" +
            "CAD3087\t\t\t\t\t\t190\tA\t169\t246\t\t\t\n" +
            "CAD3242\t\t\t\t\t\t190\tA\t173\t\t\t\t\n" +
            "CAD3428\t\t\t\t134\t186\t188\tA\t169\t246\t\t\t\n" +
            "CAD3553\t\t\t\t134\t186\t188\tA\t169\t248\t237\t\t\n" +
            "CAD3645\t\t\t\t\t182\t188\tA\t169\t248\t225\t\t\n" +
            "CAD3898\t\t\t\t\t188\t188\tA\t169\t248\t223\t\t\n" +
            "CAD4978\t\t\t\t\t188\t188\tA\t169\t248\t235\t244\t\n" +
            "8053\t\t\t\t\t184\t188\tA\t169\t246\t\t\t\n" +
            "3908\t\t\t\t\t182\t188\tA\t163\t\t\t\t\n" +
            "L2E05\t\t\t\t\t188\t188\tA\t169\t250\t\t\t\n" +
            "CAD5914\t\t\t\t\t188\t188\tA\t169\t258\t\t\t\n" +
            "NW1\t169\t177\t237\t0\t186\t188\tA\t161\t\t\t\t\n" +
            "NW3\t\t\t\t138\t186\t188\tA\t161\t\t\t\t\n" +
            "Positions\t???:???\t???:???\t???:???\t???:???\t???:???\t???:???\t55044016\t???:???\t???:???\t???:???\t???:???\t???:???";

    String expectedResultsFF = ".\tLeft3\tLeft2\tLeft1\tTarget\tRight1\tRight2\tRight3\tRight4\tRight5\n" +
            "Samples\\Markers\tD1S197\tD1S2661\tD1S417\trs28942111\tD1S200\tD1S2742\tD1S220\tD1S473\tD1S390\n" +
            "HC2\t132/136\t186\t188\tA\t169\t248\t235\t236/238\t\n" +
            "HC806\t\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "HC92\t\t182\t188\tA\t169\t248\t223\t\t\n" +
            "HC748\t\t182\t188\tA\t169\t260\t\t\t\n" +
            "HC2062\t\t0\t188\tA\t169\t248\t235\t240\t208\n" +
            "Lyon1\t132\t186\t188\tA\t169\t246\t\t\t\n" +
            "Lyon2\t142\t186\t188\tA\t169\t246\t\t\t\n" +
            "S Afr\t\t182\t188\tA\t167\t\t\t\t\n" +
            "CAD3077\t136\t186\t188\tA\t175\t\t\t\t\n" +
            "CAD3087\t\t\t190\tA\t169\t246\t\t\t\n" +
            "CAD3242\t\t\t190\tA\t173\t\t\t\t\n" +
            "CAD3428\t134\t186\t188\tA\t169\t246\t\t\t\n" +
            "CAD3553\t134\t186\t188\tA\t169\t248\t237\t\t\n" +
            "CAD3645\t\t182\t188\tA\t169\t248\t225\t\t\n" +
            "CAD3898\t\t188\t188\tA\t169\t248\t223\t\t\n" +
            "CAD4978\t\t188\t188\tA\t169\t248\t235\t244\t\n" +
            "8053\t\t184\t188\tA\t169\t246\t\t\t\n" +
            "3908\t\t182\t188\tA\t163\t\t\t\t\n" +
            "L2E05\t\t188\t188\tA\t169\t250\t\t\t\n" +
            "CAD5914\t\t188\t188\tA\t169\t258\t\t\t\n" +
            "NW1\t0\t186\t188\tA\t161\t\t\t\t\n" +
            "NW3\t138\t186\t188\tA\t161\t\t\t\t\n" +
            "Positions\t???:???\t???:???\t???:???\t55044016\t???:???\t???:???\t???:???\t???:???\t???:???";

    int col = 7; //7th marker, 1-based
    String position = "55044016";

    boolean success = true;
    success &= testPhasing(inputData, expectedResultsTF, col, position, true, false);
    success &= testPhasing(inputData, expectedResultsTT, col, position, true, true);
    success &= testPhasing(inputData, expectedResultsFT, col, position, false, true);
    success &= testPhasing(inputData, expectedResultsFF, col, position, false, false);

    if(success)
      System.err.println("[SUCCESS] Everything went well");
    else
      System.err.println("[FAILURE] At least one test failed");
  }

  private static boolean testPhasing(String inputData, String expectedResults, int col, String position, boolean stopOnExAequo, boolean ignoreMissing) throws IOException {
    File tmpDir = new File(System.getProperty("java.io.tmpdir"));
    File inputFile = File.createTempFile("test", ".estiage", tmpDir);
    inputFile.deleteOnExit();
    File outputFile = new File(inputFile+".phased");
    outputFile.deleteOnExit();
    System.err.println("Input File : "+inputFile.getAbsolutePath());
    PrintWriter out = new PrintWriter(new FileWriter(inputFile));
    out.print(inputData);
    out.flush();
    out.close();

    try {
      Main.phase(inputFile.toString(), outputFile.toString(), col, position, stopOnExAequo, ignoreMissing);
    } catch(Exception e){
      e.printStackTrace();
    }

    UniversalReader in = new UniversalReader(outputFile.toString());
    String line;
    ArrayList<String> lines = new ArrayList<>();
    while((line = in.readLine()) != null) {
      System.out.println(line);
      lines.add(line);
    }
    in.close();
    String actualOutput = String.join("\n", lines);
    boolean success = expectedResults.equals(actualOutput);
    if(success)
      System.err.println("SUCCESS : actual output matches expected output");
    else {
      String[] expect = expectedResults.split("\n",- 1);
      System.err.println("FAILURE : actual output does not match expect output");
      System.err.println("Number of actual/expected lines : "+lines.size()+"/"+expect.length);
      for(int i = 0; i < Math.max(lines.size(), expect.length); i++){
        String thisAc = "";
        String thisEx = "";
        if(i < lines.size())
          thisAc = lines.get(i);
        if(i < expect.length)
          thisEx = expect[i];
        if(!thisEx.equals(thisAc)){
          System.err.println("Line "+(i+1));
          System.err.println("Expected |"+thisEx);
          System.err.println("Actual   |"+thisAc);
        }
      }
    }
    return success;
  }

  private static void testF508Del() throws Exception {
    final String inputdata = "24 16 12\n" +
            "2.771067021944659E-4 6.265502119854636E-4 7.099501314594511E-4 7.125775084720532E-4 8.908490785231561E-4 9.15740123327412E-4 0.0021676843556302527 0.002681018607915704 0.0031482097226622134 0.003474526581724931 0.004693223819911853 0.005813958849841483 0.0063223414981944115 0.007379135971942818 0.008982089166839057 0.008983564884664312\n" +
            "0.0823529 0.502353 0.766471 0.74 0.0176471 0.9882353 0.681765 0.9882353 0.9564706 0.00588235 0.9723529 0.534706 0.628235 0.418235 0.9835294 0.0\n" +
            "1.5765061846132378E-4 2.8372539595914034E-4 2.838634824378612E-4 3.592323267671148E-4 5.837335236468365E-4 7.621413500345191E-4 8.272063539536219E-4 0.0018748766000957983 0.0018886781913002788 0.00246193731359268 0.004938675059907856 0.019097818871431478\n" +
            "0.0211765 0.849412 0.99705882 0.9876471 0.611176 0.616471 0.388235 0.269412 0.577059 0.375882 0.0135294 0.0\n" +
            "0.001 0\n" +
            "1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -1\n" +
            "1 1 1 1 1 1 1 1 1 1 1 -1\n" +
            "5 5 2 2\n" +
            "10 6 2 2\n" +
            "5 7 2 2\n" +
            "13 11 2 2\n" +
            "2 7 2 2\n" +
            "4 7 2 2\n" +
            "2 3 2 2\n" +
            "2 3 2 2\n" +
            "16 9 1 2\n" +
            "11 12 2 1\n" +
            "13 1 2 2\n" +
            "7 7 2 2\n" +
            "10 7 2 2\n" +
            "1 3 2 2\n" +
            "6 8 2 2\n" +
            "16 12 1 1\n" +
            "14 9 2 2\n" +
            "2 3 2 2\n" +
            "10 6 2 2\n" +
            "12 9 2 2\n" +
            "11 10 2 2\n" +
            "9 7 2 2\n" +
            "8 4 2 2\n" +
            "3 3 2 2";
    final String expectedResults = "n = 149, nend = 294, ninf = 109, nsup = 209, likelihood = -114.29372603629017";
    //write to tmp
    File tmpFile = Files.createTempFile("test", ".estiage").getFileName().toFile();
    PrintWriter out = new PrintWriter(new FileWriter(tmpFile));
    out.print(inputdata);
    out.flush();
    out.close();
    System.out.println("Expected Results:");
    System.out.println(expectedResults);
    System.out.println("Actual Results:");
    Main.run(tmpFile.toString());

    Files.delete(tmpFile.toPath());
  }
}
