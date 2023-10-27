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
  public static final String KEY_PHASE = "phase";
  public static final String KEY_VCF2COMPLETE = "vcf2complete";
  public static final String KEY_VCF2PREINPUT = "vcf2preinput";
  public static final String KEY_VCF2INPUT = "vcf2input";
  public static final String KEY_RAW2COMPLETE = "raw2complete";
  public static final String KEY_RAW2PREINPUT = "raw2preinput";
  public static final String KEY_RAW2INPUT = "raw2input";
  public static final String KEY_COMPLETE2PREINPUT = "complete2preinput";
  public static final String KEY_COMPLETE2INPUT = "complete2input";
  public static final String KEY_PREINPUT2INPUT = "preinput2input";
  public static final String KEY_RUN = "run";
  public static final String KEY_RATE = "rate";
  public static final String KEY_NO_COLOR = "--nocolor";


  public static final String EXT_PHASED = ".phased";
  public static final String EXT_RAW = ".estiraw";
  public static final String EXT_FULL = ".estifull";
  public static final String EXT_PREINPUT = ".preinput";
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

    String vcf, raw, complete, preinput, input, chrPosAllele, gnomad, hapmap, mutationModel, mutationRate;
    VCFFile.Mode mode;

    switch(args[0].toLowerCase()){
      case KEY_RATE:
        if(args.length < 4)
          usagerate(true);
        String hapmapfile = args[1];
        String p1 = args[2];
        String p2 = args[3];
        rate(hapmapfile, p1, p2);
        break;
      case KEY_PHASE:
        if(args.length < 3)
          usagephase(true);
        String inputfile = args[1];
        String col = args[2];
        String position = args[3];
        int colnum = -1;
        try{
          colnum = Integer.parseInt(col);
        } catch(NumberFormatException e){
          usagephase(true);
        }
        phase(inputfile, colnum, position);
        break;
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
      case KEY_COMPLETE2PREINPUT:
        if(args.length < 5)
          usagecomplete2preinput(true);
        complete = args[1];
        preinput = args[2];
        mutationModel = args[3];
        mutationRate = args[4];
        complete2preinput(complete, preinput, mutationModel, mutationRate);
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
      case KEY_RAW2PREINPUT:
        if(args.length < 7)
          usageraw2preinput(true);
        raw = args[1];
        preinput = args[2];
        gnomad = args[3];
        hapmap = args[4];
        mutationModel = args[5];
        mutationRate = args[6];
        raw2preinput(raw, preinput, gnomad, hapmap, mutationModel, mutationRate);
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
      case KEY_VCF2PREINPUT:
        if(args.length < 9)
          usagevcf2preinput(true);
        vcf = args[1];
        preinput = args[2];
        chrPosAllele = args[3];
        mode = VCFFile.Mode.valueOf(args[4]);
        gnomad = args[5];
        hapmap = args[6];
        mutationModel = args[7];
        mutationRate = args[8];
        vcf2preinput(vcf, preinput, chrPosAllele, mode, gnomad, hapmap, mutationModel, mutationRate);
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
      case KEY_PREINPUT2INPUT:
        if(args.length < 3)
          usagepreinput2input(true);
        preinput = args[1];
        input = args[2];
        preinput2input(preinput, input);
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
    usagephase(false);
    usagevcf2raw(false);
    usageraw2complete(false);
    usagevcf2complete(false);
    usagecomplete2input(false);
    usageraw2input(false);
    usagevcf2input(false);
    usagecomplete2preinput(false);
    usageraw2preinput(false);
    usagevcf2preinput(false);
    usagepreinput2input(false);
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

  private static void usagerate(boolean printPrefix){
    printUsage(printPrefix, KEY_RATE, "HapMapFilename", "Position1", "Position2");
  }

  private static void usagephase(boolean printPrefix){
    printUsage(printPrefix, KEY_PHASE, INPUT, INPUT+EXT_PHASED, "Column", "TargetPosition");
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

  private static void usagecomplete2preinput(boolean printPrefix){
    printUsage(printPrefix, KEY_COMPLETE2PREINPUT, INPUT+EXT_FULL, OUTPUT+EXT_PREINPUT, MODELS, RATE);
  }

  private static void usageraw2preinput(boolean printPrefix){
    printUsage(printPrefix, KEY_RAW2INPUT, INPUT+EXT_RAW, OUTPUT+EXT_PREINPUT, GNOMAD, HAPMAP, MODELS, RATE);
  }

  private static void usagevcf2preinput(boolean printPrefix){
    printUsage(printPrefix, KEY_VCF2INPUT, INPUT+EXT_VCF, OUTPUT+EXT_PREINPUT, CHROMPOSALLELE, VCFMODE, GNOMAD, HAPMAP, MODELS, RATE);
  }

  private static void usagepreinput2input(boolean printPrefix){
    printUsage(printPrefix, KEY_PREINPUT2INPUT, INPUT+EXT_PREINPUT, OUTPUT+EXT_INPUT);
  }

  private static void usagerun(boolean printPrefix){
    printUsage(printPrefix, KEY_RUN, INPUT+EXT_INPUT);
  }

  private static void vcf2raw(String vcf, String raw, String chrPosAllele, VCFFile.Mode mode) throws InterruptedException, EstiageFormatException, IOException {
    VCFFile vcfFile = new VCFFile(vcf, mode);
    vcfFile.setVariant(chrPosAllele);
    vcfFile.exportAsRaw(raw);
  }

  private static void rate(String hapMapFilename, String pos1, String pos2) throws IOException, EstiageFormatException {
    String hapmapfile = hapMapFilename;
    int p1 = Integer.parseInt(pos1);
    int p2 = Integer.parseInt(pos2);
    int distance = 1 + p2 - p1;
    double mb = distance*0.000001;
    HapMap hapMap = new HapMap(hapmapfile, p1, p2);
    System.out.println("Measuring recombination fraction between ["+p1+"] and ["+p2+"] from file : "+hapmapfile);

    double rate = hapMap.getRate(p1, p2);
    double cM = MathLib.getCentiMorgan(rate, mb);
    double d = cM / 100; //from cMorgans to Morgans
    double recombinationFraction = MathLib.kosambiTheta(d);
    System.out.println("distance (b) : "+distance);
    System.out.println("distance (Mb) : "+mb);
    System.out.println("Rate : "+rate);
    System.out.println("cM : "+cM);
    System.out.println("Morgans : "+d);
    System.out.println("theta : "+recombinationFraction);
  }

  private static void phase(String inputfile, int col, String position) throws InterruptedException, EstiageFormatException, IOException {
    Unphased unphased = new Unphased(inputfile, col, position);
    unphased.phase();
    unphased.export(inputfile + EXT_PHASED);
    /*
    VCFFile vcfFile = new VCFFile(vcf, mode);
    vcfFile.setVariant(chrPosAllele);
    vcfFile.exportAsRaw(raw);
    */
  }

  private static void raw2complete(String raw, String complete, String gnomad, String hapmap) throws IOException, EstiageFormatException {
    TSVFile rawfile = new TSVFile(raw, TSVFile.Type.RAW);
    rawfile.printSummary();
    //Getting frequencies
    rawfile.generateFrequenciesAndRecombinationFractions(gnomad, hapmap);
    rawfile.export(complete);
  }

  private static void complete2preinput(String complete, String preinput, String mutationModel, String mutationRate) throws IOException, EstiageFormatException {
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
    estiageInput.printToFile(preinput);
  }

  private static void complete2input(String complete, String input, String mutationModel, String mutationRate) throws IOException, EstiageFormatException {
    String preinput = input.replace(EXT_INPUT, EXT_PREINPUT);
    complete2preinput(complete, preinput, mutationModel, mutationRate);
    preinput2input(preinput, input);
  }

  private static void preinput2input(String preinput, String input) throws IOException, EstiageFormatException {
    InputFile iFile = new InputFile(preinput);
    iFile.fromPreinput2Input();
    iFile.printToFile(input);
  }

  private static void vcf2complete(String vcf, String complete, String chrPosAllele, VCFFile.Mode mode, String gnomad, String hapmap) throws IOException, EstiageFormatException, InterruptedException {
    String raw = vcf+EXT_RAW;
    vcf2raw(vcf, raw, chrPosAllele, mode);
    raw2complete(raw, complete, gnomad, hapmap);
  }

  private static void raw2preinput(String raw, String preinput, String gnomad, String hapmap, String mutationModel, String mutationRate) throws IOException, EstiageFormatException {
    String complete = raw+EXT_FULL;
    raw2complete(raw, complete, gnomad, hapmap);
    complete2preinput(complete, preinput, mutationModel, mutationRate);
  }

  public static void raw2input(String raw, String input, String gnomad, String hapmap, String mutationModel, String mutationRate) throws IOException, EstiageFormatException {
    String preinput = input.replace(EXT_INPUT, EXT_PREINPUT);
    raw2preinput(raw, preinput, gnomad, hapmap, mutationModel, mutationRate);
    preinput2input(preinput, input);
  }

  private static void vcf2preinput(String vcf, String preinput, String chrPosAllele, VCFFile.Mode mode, String gnomad, String hapmap, String mutationModel, String mutationRate) throws IOException, EstiageFormatException, InterruptedException {
    String raw = vcf+EXT_RAW;
    vcf2raw(vcf, raw, chrPosAllele, mode);
    raw2preinput(raw, preinput, gnomad, hapmap, mutationModel, mutationRate);
  }

  private static void vcf2input(String vcf, String input, String chrPosAllele, VCFFile.Mode mode, String gnomad, String hapmap, String mutationModel, String mutationRate) throws IOException, EstiageFormatException, InterruptedException {
    String preinput = input.replace(EXT_INPUT, EXT_PREINPUT);
    vcf2preinput(vcf, preinput, chrPosAllele, mode, gnomad, hapmap, mutationModel, mutationRate);
    preinput2input(preinput, input);
  }

  public static void run(String filename) throws IOException, EstiageCTranslation.EstiageException {
    EstiageCTranslation e = new EstiageCTranslation(filename);
    e.run();
  }
}
