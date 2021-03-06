package fr.inserm.u1078.estiage;

import fr.inserm.u1078.tludwig.maok.tools.Message;

import java.io.IOException;

/**
 * Main Class for Estiage
 *
 * @author Thomas E. Ludwig (INSERM - U1078)
 * Started on             2021-03-17
 * Checked for release on XXXX-XX-XX
 * Unit Test defined on   XXXX-XX-XX
 */
public class Main {

  public static final String KEY_VCF2RAW = "vcf2raw";
  public static final String KEY_VCF2COMPLETE = "vcf2complete";
  public static final String KEY_VCF2INPUT = "vcf2input";
  public static final String KEY_RAW2COMPLETE = "raw2complete";
  public static final String KEY_RAW2INPUT = "raw2input";
  public static final String KEY_COMPLETE2INPUT = "complete2input";
  public static final String KEY_RUN = "run";
  public static final String KEY_NO_COLOR = "--nocolor";

  public static final String EXT_RAW = ".estiraw";
  public static final String EXT_FULL = ".estifull";
  public static final String EXT_INPUT = ".estinput";
  public static final String EXT_VCF = ".vcf(.gz)";

  public static final String MODELS = "mutationModel[0:normal|1:stepwise]";
  public static final String RATE = "mutationRate";
  public static final String GNOMAD = "Gnomad"+EXT_VCF;
  public static final String HAPMAP = "hapmap.txt";
  public static final String CHROMPOSALLELE = "chr:pos:allele(homoAltAllele)";
  public static final String INPUT = "input";
  public static final String OUTPUT = "output";
  public static final String ESTIAGE = "EstiAge";
  public static final String VCFMODE = "mode["+VCFFile.Mode.IGNORE+"|"+VCFFile.Mode.HETEROZYGOUS+"|"+VCFFile.Mode.HOMOZYGOUS+"]";

  public static void main(String[]args) throws IOException, EstiageFormatException, EstiageCTranslation.EstiageException, InterruptedException {
    Message.setDebugActive(true);
    for(String arg : args)
      if(KEY_NO_COLOR.equalsIgnoreCase(arg))
        Message.setWithColor(false);

    StringBuilder message = new StringBuilder("Running with arguments:");
    for(String arg : args)
      message.append(" [").append(arg).append("]");
    message.append(".");
    Message.info(message.toString());
    if(args.length < 1)
      usage();

    String vcf, raw, complete, input, chrPosAllele, gnomad, hapmap, mutationModel, mutationRate;
    VCFFile.Mode mode;

    switch(args[0].toLowerCase()){
      case KEY_VCF2RAW:
        if(args.length < 5)
          usagevcf2raw(true);
        vcf = args[1];
        raw = args[2];
        chrPosAllele = args[3];
        mode = VCFFile.Mode.valueOf(args[4].toUpperCase());
        vcf2raw(vcf, raw, chrPosAllele, mode);
        break;
      case KEY_RAW2COMPLETE:
        if(args.length < 5)
          usageraw2complete(true);
        raw = args[1];
        complete = args[2];
        gnomad = args[3];
        hapmap = args[4];
        raw2complete(raw, complete, gnomad, hapmap);
        break;
      case KEY_VCF2COMPLETE:
        if(args.length < 7)
          usagevcf2complete(true);
        vcf = args[1];
        complete = args[2];
        chrPosAllele = args[3];
        mode = VCFFile.Mode.valueOf(args[4]);
        gnomad = args[5];
        hapmap = args[6];
        vcf2complete(vcf, complete, chrPosAllele, mode, gnomad, hapmap);
        break;
      case KEY_COMPLETE2INPUT:
        if(args.length < 5)
          usagecomplete2input(true);
        complete = args[1];
        input = args[2];
        mutationModel = args[3];
        mutationRate = args[4];
        complete2input(complete, input, mutationModel, mutationRate);
        break;
      case KEY_RAW2INPUT:
        if(args.length < 7)
          usageraw2input(true);
        raw = args[1];
        input = args[2];
        gnomad = args[3];
        hapmap = args[4];
        mutationModel = args[5];
        mutationRate = args[6];
        raw2input(raw, input, gnomad, hapmap, mutationModel, mutationRate);
        break;
      case KEY_VCF2INPUT:
        if(args.length < 9)
          usagevcf2input(true);
        vcf = args[1];
        input = args[2];
        chrPosAllele = args[3];
        mode = VCFFile.Mode.valueOf(args[4]);
        gnomad = args[5];
        hapmap = args[6];
        mutationModel = args[7];
        mutationRate = args[8];
        vcf2input(vcf, input, chrPosAllele, mode, gnomad, hapmap, mutationModel, mutationRate);
        break;
      case KEY_RUN:
        if(args.length < 2)
          usagerun(true);
        String filename = args[1];
        run(filename);
        break;
      default :
        Message.error("Unknown Option ["+args[0]+"]");
        usage();
    }
  }

  private static void usage(){
    System.err.println(ESTIAGE+"\nUsage :");
    usagevcf2raw(false);
    usageraw2complete(false);
    usagevcf2complete(false);
    usagecomplete2input(false);
    usageraw2input(false);
    usagevcf2input(false);
    usagerun(false);

    System.exit(1);
  }

  private static void printUsage(boolean printPrefix, String... args){
    if(printPrefix)
      System.err.println(ESTIAGE+"\nUsage :");
    System.err.println("\t"+String.join(" ",args));
    if(printPrefix)
      System.exit(1);
  }

  private static void usagevcf2raw(boolean printPrefix){
    printUsage(printPrefix, KEY_VCF2RAW, INPUT+EXT_VCF, OUTPUT+EXT_RAW, CHROMPOSALLELE, VCFMODE);
  }

  private static void usageraw2complete(boolean printPrefix){
    printUsage(printPrefix, KEY_RAW2COMPLETE, INPUT+EXT_RAW, OUTPUT+EXT_FULL, GNOMAD, HAPMAP);
  }

  private static void usagevcf2complete(boolean printPrefix){
    printUsage(printPrefix, KEY_VCF2COMPLETE, INPUT+EXT_VCF, OUTPUT+EXT_FULL, CHROMPOSALLELE, VCFMODE, GNOMAD, HAPMAP);
  }

  private static void usagecomplete2input(boolean printPrefix){
    printUsage(printPrefix, KEY_COMPLETE2INPUT, INPUT+EXT_FULL, OUTPUT+EXT_INPUT, MODELS, RATE);
  }

  private static void usageraw2input(boolean printPrefix){
    printUsage(printPrefix, KEY_RAW2INPUT, INPUT+EXT_RAW, OUTPUT+EXT_INPUT, GNOMAD, HAPMAP, MODELS, RATE);
  }

  private static void usagevcf2input(boolean printPrefix){
    printUsage(printPrefix, KEY_VCF2INPUT, INPUT+EXT_VCF, OUTPUT+EXT_INPUT, CHROMPOSALLELE, VCFMODE, GNOMAD, HAPMAP, MODELS, RATE);
  }

  private static void usagerun(boolean printPrefix){
    printUsage(printPrefix, KEY_RUN, INPUT+EXT_INPUT);
  }

  private static void vcf2raw(String vcf, String raw, String chrPosAllele, VCFFile.Mode mode) throws InterruptedException, EstiageFormatException, IOException {
    VCFFile vcfFile = new VCFFile(vcf, mode);
    vcfFile.setVariant(chrPosAllele);
    vcfFile.exportAsRaw(raw);
  }

  private static void raw2complete(String raw, String complete, String gnomad, String hapmap) throws IOException, EstiageFormatException {
    TSVFile rawfile = new TSVFile(raw, TSVFile.Type.RAW);
    rawfile.printSummary();
    //Getting frequencies
    rawfile.generateFrequenciesAndRecombinationFractions(gnomad, hapmap);
    rawfile.export(complete);
  }

  private static void complete2input(String complete,  String input, String mutationModel, String mutationRate) throws IOException, EstiageFormatException {
    int model = -1;
    double rate = 0;
    try {
      model = Integer.parseInt(mutationModel);
    } catch(NumberFormatException ignore) {
      //Nothing
    }
    if(model < 0 || model > 1) {
      System.err.println("Unexpected Mutation Model ["+mutationModel+"]. Should be 0 (normal) or 1 (stepwise)");
      System.exit(1);
    }
    try{
      rate = Double.parseDouble(mutationRate);
    } catch(NumberFormatException e){
      System.err.println("Could not parse mutation rate ["+mutationRate+"]");
      System.exit(1);
    }
    Message.info("Opening complete file");
    TSVFile completeFile = new TSVFile(complete, TSVFile.Type.COMPLETE);
    Message.info("Building Estiage input file");
    InputFile estiageInput = new InputFile(completeFile, model, rate);
    Message.info("Writing filename");
    estiageInput.printToFile(input);
  }

  private static void vcf2complete(String vcf, String complete, String chrPosAllele, VCFFile.Mode mode, String gnomad, String hapmap) throws IOException, EstiageFormatException, InterruptedException {
    String raw = vcf+EXT_RAW;
    vcf2raw(vcf, raw, chrPosAllele, mode);
    raw2complete(raw, complete, gnomad, hapmap);
  }

  private static void raw2input(String raw, String input, String gnomad, String hapmap, String mutationModel, String mutationRate) throws IOException, EstiageFormatException {
    String complete = raw+EXT_FULL;
    raw2complete(raw, complete, gnomad, hapmap);
    complete2input(complete, input, mutationModel, mutationRate);
  }

  private static void vcf2input(String vcf, String input, String chrPosAllele, VCFFile.Mode mode, String gnomad, String hapmap, String mutationModel, String mutationRate) throws IOException, EstiageFormatException, InterruptedException {
    String raw = vcf+EXT_RAW;
    vcf2raw(vcf, raw, chrPosAllele, mode);
    raw2input(raw, input, gnomad, hapmap, mutationModel, mutationRate);
  }

  private static void run(String filename) throws IOException, EstiageCTranslation.EstiageException {
    EstiageCTranslation e = new EstiageCTranslation(filename);
    e.run();
  }
}
