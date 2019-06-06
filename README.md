# Maize by Michelangelo
Pipeline for designing whole-chromosome oligo-FISH paints. Inspired by Jiming Jiang and Jim Birchler's [chromosome-specific FISH library in maize](https://doi.org/10.1073/pnas.1813957116) (codename: "Maize by Monet").

Developed by Lisa Malins under direction from Kirk Amundson in the [Comai Lab](http://comailab.genomecenter.ucdavis.edu/index.php/Main_Page) at the [UC Davis Genome Center](https://genomecenter.ucdavis.edu/).

## Steps
1. Divide a reference genome into oligos: for example all 45-mers separated by a step size of 3.
2. Map these oligos back to the reference genome to see where they would hybridize and how well.
3. Filter out oligos that map back to more than one location.
4. Filter out oligos that would make poor probes because they would not hybridize well or would form hairpins.
5. Calculate k-mer scores for each oligo kept in the pipeline.

## References
Albert PS, Zhang T, Semrau K, Rouillard JM, Kao YH, Wang CJ, Danilova TV, Jiang J, Birchler JA. Whole-chromosome paints in maize reveal rearrangements, nuclear domains, and chromosomal relationships. *Proceedings of the National Academy of Sciences*. 2019 Jan 29;116(5):1679-85.
