package org.monarchinitiative;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.monarchinitiative.phenol.base.PhenolRuntimeException;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.phenopackets.schema.v1.core.Disease;
import org.phenopackets.schema.v1.core.Gene;
import org.phenopackets.schema.v1.core.OntologyClass;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

import static org.monarchinitiative.phenol.formats.hpo.HpoModeOfInheritanceTermIds.*;

public class CountPhenoP {
    private static final Logger logger = LoggerFactory.getLogger(CountPhenoP.class);
    private String phenopacketDirectoryPath;
    private List<File> phenopacketFiles;
    private Ontology ontology;
    private Map<TermId, HpoDisease> diseaseMap;

    private Map<TermId, Integer> disease2count = new HashMap<>();

    private Map<TermId, Integer> hpo2count = new HashMap<>();

    private Set<String> genes = new HashSet<>();

    private DescriptiveStatistics termsPerPhenopacket = new DescriptiveStatistics();
    private DescriptiveStatistics negatedTermsPerPhenopacket = new DescriptiveStatistics();

    double medianCountPerDisease;
    double maxCountPerDiseases;
    int n_recessive =0;
    int n_dominant = 0;
    int n_xchromosomal = 0;
    int n_heterogeneous = 0;
    int n_somatic = 0;
    int n_sporadic = 0;
    int n_somatic_mosaic = 0;


    public static void main(String []args) throws IOException {
        String hpPath = "/home/robinp/IdeaProjects/LIRICAL/data/hp.obo";
        String ppacketDir = "/home/robinp/Desktop/ppacket";
        String phenotypeAnnotationPath = "/home/robinp/IdeaProjects/LIRICAL/data/phenotype.hpoa";
        CountPhenoP cpp = new CountPhenoP(hpPath,ppacketDir, phenotypeAnnotationPath);
        BufferedWriter writer = new BufferedWriter(new FileWriter("phenopacketstats.tex"));
        cpp.writeLongTable(writer);
        cpp.getStats(writer);
        cpp.printStats(writer);

        writer.write("\\end{longtable}\n");
        writer.close();
    }



    public CountPhenoP(String hpoPath, String ppacketDirPath, String phenotypeAnnotationPath) {
        this.ontology = OntologyLoader.loadOntology(new File(hpoPath));
        diseaseMap = HpoDiseaseAnnotationParser.loadDiseaseMap(phenotypeAnnotationPath,ontology);
        this.phenopacketDirectoryPath = ppacketDirPath;
        getListOfPhenopacketFiles();
    }


    private void recordPhenotypes(List<TermId> ids, List<TermId> negated) {
        termsPerPhenopacket.addValue(ids.size());
        negatedTermsPerPhenopacket.addValue(negated.size());
        for (TermId tid : ids) {
            this.hpo2count.putIfAbsent(tid,0);
            this.hpo2count.put(tid,this.hpo2count.get(tid));
        }
        for (TermId tid : negated) {
            this.hpo2count.putIfAbsent(tid,0);
            this.hpo2count.put(tid,1+this.hpo2count.get(tid));
        }
    }

    private void recordDiagnosis(Disease d) {
        OntologyClass oc = d.getTerm();
        TermId did = TermId.of(oc.getId());
        this.disease2count.putIfAbsent(did,0);
        this.disease2count.put(did,1 + disease2count.get(did));
        HpoDisease hpod = diseaseMap.get(did);
        if(hpod == null) {
            System.out.println("Could not retrieve data for " + d.toString());
            return;
        }
        List<TermId> inheritance = hpod.getModesOfInheritance();
        if (inheritance == null) {
            System.out.println("Could not find inheritance term for " + hpod.getName());
            return;
        }
        for (TermId id : inheritance) {
            if (id.equals(AUTOSOMAL_RECESSIVE)) {
                n_recessive++;
            } else if (id.equals(AUTOSOMAL_DOMINANT) || id.equals(CONTIGUOUS_GENE_SYNDROME_AUTOSOMAL_DOMINANT)) {
                n_dominant++;
            } else if (id.equals(X_LINKED) || id.equals(X_LINKED_RECESSIVE) || id.equals(X_LINKED_DOMINANT)) {
                n_xchromosomal++;
            } else if (id.equals(HETEROGENEOUS)) {
                n_heterogeneous++;
            } else if (id.equals(SOMATIC_MUTATION)) {
                n_somatic++;
            } else if (id.equals(SPORADIC)) {
                n_sporadic++;
            }  else if (id.equals(SOMATIC_MOSAICISM)) {
                n_somatic_mosaic++;
            }  else{
                System.out.println("Could not identify id " + id.getValue());
                System.exit(0);
            }
        }

    }


    public void printStats(Writer writer) {
        System.out.printf("Number of diseases: %d (median %f, max %f)\n", this.disease2count.size(), medianCountPerDisease,maxCountPerDiseases);
        System.out.printf("Autosomal recessive: %d\n", this.n_recessive);
        System.out.printf("Autosomal dominant: %d\n", this.n_dominant);
        System.out.printf("X chromosomal: %d\n", this.n_xchromosomal);
        System.out.printf("heterogeneous: %d\n", this.n_heterogeneous);
        System.out.printf("somatic: %d\n", this.n_somatic);
        System.out.printf("somatic mosaic: %d\n", this.n_somatic_mosaic);
        System.out.printf("sporadic: %d\n", this.n_somatic);
        System.out.printf("Number of genes: %d\n", genes.size());
        System.out.printf("Total number of HPO terms used in phenopackets: %d\n", this.hpo2count.size());
        int mn = 0;
        for (Integer i : this.hpo2count.values()) {
            mn += i;
        }
        System.out.printf("Mean number of times each HPO term was used: %.2f\n", (double)mn/this.hpo2count.size());
        int medianNegated = (int)negatedTermsPerPhenopacket.getPercentile(50.0);
        double meanNegated = negatedTermsPerPhenopacket.getMean();
        int maxNegated = (int)negatedTermsPerPhenopacket.getMax();
        System.out.printf("Mean negated %f, median %d max %d\n",meanNegated, medianNegated,maxNegated);
        int median = (int)termsPerPhenopacket.getPercentile(50.0);
        double mean = termsPerPhenopacket.getMean();
        int max = (int)negatedTermsPerPhenopacket.getMax();
        System.out.printf("Mean  %f, median %d max %d\n",mean, median,max);
    }


    public void writeLongTable(Writer writer) throws IOException {
        writer.write("\\begin{longtable}{|l|r|r|r|l|}\n" +
                "\\caption{Phenopackets analyzed in this work}\\\\\n" +
                "\\hline\n" +
                "\\textbf{Disease} & \\textbf{Gene} & \\textbf{Proband} & \\textbf{n. HPO terms}& \\textbf{Publication} \\\\\n" +
                "\\hline\n" +
                "\\endfirsthead\n" +
                "\\multicolumn{5}{c}%\n" +
                "{\\tablename\\ \\thetable\\ -- \\textit{Continued from previous page}} \\\\\n" +
                "\\hline\n" +
                "\\textbf{Disease} & \\textbf{Gene} & \\textbf{Proband} & \\textbf{n. HPO terms}& \\textbf{Publication} \\\\\n" +
                "\\hline\n" +
                "\\endhead\n" +
                "\\hline \\multicolumn{5}{r}{\\textit{Continued on next page}} \\\\\n" +
                "\\endfoot\n" +
                "\\hline\n" +
                "\\endlastfoot");
    }

    static String convert(String str)
    {

        // Create a char array of given String
        char ch[] = str.toCharArray();
        for (int i = 0; i < str.length(); i++) {

            // If first character of a word is found
            if (i == 0 && ch[i] != ' ' ||
                    ch[i] != ' ' && (ch[i - 1] == ' ' || ch[i-1] == '-')) {

                // If it is in lower-case
                if (ch[i] >= 'a' && ch[i] <= 'z') {

                    // Convert into Upper-case
                    ch[i] = (char)(ch[i] - 'a' + 'A');
                }
            }

            // If apart from first character
            // Any one is in Upper-case
            else if (ch[i] >= 'A' && ch[i] <= 'Z')

                // Convert into Lower-Case
                ch[i] = (char)(ch[i] + 'a' - 'A');
        }

        // Convert the char array to equivalent String
        String st = new String(ch);
        return st;
    }

    // disease, gene, proband, hpoterms, pub.
    public void getStats(Writer writer)   throws IOException {
        for (java.io.File file: this.phenopacketFiles) {
            if (! file.exists()) {
                throw new RuntimeException("Could not find phenopacket file at " + file.getAbsolutePath());
            }
            PhenopacketImporter importer = PhenopacketImporter.fromJson(file.getAbsolutePath(), this.ontology);
            Disease disease = importer.getDiagnosis();
            recordDiagnosis(disease);
            List<TermId> ids = importer.getHpoTerms();
            List<TermId> negated = importer.getNegatedHpoTerms();
            recordPhenotypes(ids,negated);
            Gene g = importer.getGene();
            genes.add(g.getId());
            String diseaseName = disease.getTerm().getLabel();
            int i = diseaseName.indexOf(";");
            if (i>0) {
                diseaseName = diseaseName.substring(0,i);
            }
            diseaseName = convert(diseaseName);
            diseaseName = diseaseName.replace("Syndrome", "syndrome");
            writer.write(String.format("%s & %s & %s & %d & %s\\\\ \n",diseaseName,
                    g.getSymbol(), importer.getSamplename(),(ids.size() + negated.size()),
                    importer.getPMID()));
        }

        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (TermId tid : this.disease2count.keySet()) {
            int c = disease2count.get(tid);
            stats.addValue(c);
        }
        this.medianCountPerDisease = stats.getPercentile(0.5);
        this.maxCountPerDiseases = stats.getMax();

    }





    private void getListOfPhenopacketFiles() {
        phenopacketFiles = new ArrayList<>();
        final File folder = new File(phenopacketDirectoryPath);
        if (! folder.isDirectory()) {
            throw new PhenolRuntimeException("Could not open Phenopackets directory at "+phenopacketDirectoryPath);
        }
        int counter=0;
        for (final File fileEntry : folder.listFiles()) {
            if (fileEntry.isFile() && fileEntry.getAbsolutePath().endsWith(".json")) {
                logger.info("\tPhenopacket: \"{}\"", fileEntry.getAbsolutePath());
                //System.out.println(++counter + ") "+ fileEntry.getName());
                this.phenopacketFiles.add(fileEntry);
            }
        }
    }

}
