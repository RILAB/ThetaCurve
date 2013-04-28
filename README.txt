Thetacurve is a perl script for calculating the a multilocus maximum likelihood estimate of the per locus value of Watterson's theta diversity statistic.  The script uses the recursion equations of Hudson (1991), and then multiplies likelihoods across loci.  This method assumes a constant mutation rate across loci, and makes the (conservative) assumption of no recombination.  Included is an example file with locus data taken from Wright et al.(2003). The script is used as follows:

perl Thetacurve.pl < infile

where infile can be the name of any text file. The infile should be structured exactly like the example file, with the name of a locus, the number of sequences/alleles, the number of segregating sites, and finally the number of base pairs analyzed each on separate lines.

The outfile "outfile.txt" will give the maximum likelihood estimate of theta for each locus as well as the multilocus maximum likelihood estimate.  For the multilocus values, the script will aslo print out the relative log likelihood of each value of PER LOCUS theta (for graphing) and the 95% credible interval of the multilocus estimate.

For questions, please write rossibarra@gmail.com.


Hudson RR (1990) Gene genealogies and the coalescent process.
Oxford Surveys in Evolutionary Biology, 7, 1–45.

Wright, SI, et al. (2003) Subdivision and haplotype structure in natural populations of Arabidopsis lyrata.
Molecular Ecology, 12, 1247-1263.
