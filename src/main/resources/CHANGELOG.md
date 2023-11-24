# Changelog
## 1.0.3 (2023-11-24)
- `Utils` : **added** default location of tabix can be customized with call to `java -Dtabix=/PATH/TO/tabix`
- `EstiageCTranslation` : **changed** - refactoring and code simplification
- `VCFVariant` : **moved** from `VCF.Variant`
- `VCF.exportAsRaw()` : **update** Method='CLASSICAL' or 'LONGEST_HAPLOTYPE'`
- `VCF.exportAsRaw()` : **modified** when several alleles are shared between X alleles, only keep the one further away from the target
## 1.0.2 (2023-10-27)
- `Main` : Call to get recombination fraction between 2 positions
- `MathLib` : Externalisation of Math calls
- `Unphased` : Process a new data format for unphased data
- `InputFile` : Fix export to estiage format when mixing snp/indels
## 1.0.1 (2022-08-01)
- `Estiage Format`: **optimized** Check Odd/Even Microsat
- `Estiage Format`: **added** Conversion from pre-estiage format to estiage (200->202 to 1->2)
## 1.0.0 (2022-05-03)
- First Release
