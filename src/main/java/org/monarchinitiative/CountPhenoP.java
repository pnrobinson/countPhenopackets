package org.monarchinitiative;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.monarchinitiative.phenol.base.PhenolRuntimeException;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.phenopackets.schema.v1.core.Disease;
import org.phenopackets.schema.v1.core.OntologyClass;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.util.*;

import static org.monarchinitiative.phenol.formats.hpo.HpoModeOfInheritanceTermIds.*;

public class CountPhenoP {
    private static final Logger logger = LoggerFactory.getLogger(CountPhenoP.class);
    private String phenopacketDirectoryPath;
    private List<File> phenopacketFiles;
    private Ontology ontology;
    private Map<TermId, HpoDisease> diseaseMap;

    private Map<TermId, Integer> disease2count = new HashMap<>();

    double medianCountPerDisease;
    double maxCountPerDiseases;
    int n_recessive =0;
    int n_dominant = 0;
    int n_xchromosomal = 0;
    int n_heterogeneous = 0;
    int n_somatic = 0;
    int n_sporadic = 0;
    int n_somatic_mosaic = 0;

    public static void main(String []args) {
        String hpPath = "/home/robinp/IdeaProjects/LIRICAL/data/hp.obo";
        String ppacketDir = "/home/robinp/Desktop/ppacket";
        String phenotypeAnnotationPath = "/home/robinp/IdeaProjects/LIRICAL/data/phenotype.hpoa";
        CountPhenoP cpp = new CountPhenoP(hpPath,ppacketDir, phenotypeAnnotationPath);
        cpp.getStats();
        cpp.printStats();
    }



    public CountPhenoP(String hpoPath, String ppacketDirPath, String phenotypeAnnotationPath) {
        this.ontology = OntologyLoader.loadOntology(new File(hpoPath));
        diseaseMap = HpoDiseaseAnnotationParser.loadDiseaseMap(phenotypeAnnotationPath,ontology);
        this.phenopacketDirectoryPath = ppacketDirPath;
        getListOfPhenopacketFiles();
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


    public void printStats() {
        System.out.printf("Number of diseases: %d (median %f, max %f)\n", this.disease2count.size(), medianCountPerDisease,maxCountPerDiseases);
        System.out.printf("Autosomal recessive: %d\n", this.n_recessive);
        System.out.printf("Autosomal dominant: %d\n", this.n_dominant);
        System.out.printf("X chromosomal: %d\n", this.n_xchromosomal);
        System.out.printf("heterogeneous: %d\n", this.n_heterogeneous);
        System.out.printf("somatic: %d\n", this.n_somatic);
        System.out.printf("somatic mosaic: %d\n", this.n_somatic_mosaic);
        System.out.printf("sporadic: %d\n", this.n_somatic);
    }


    public void getStats() {
        for (java.io.File file: this.phenopacketFiles) {
            String phenopacketAbsolutePath = file.getAbsolutePath();
            //File f = new File(file);
            if (! file.exists()) {
                throw new RuntimeException("Could not find phenopacket file at " + file.getAbsolutePath());
            }
            PhenopacketImporter importer = PhenopacketImporter.fromJson(file.getAbsolutePath(), this.ontology);
            Disease disease = importer.getDiagnosis();
            recordDiagnosis(disease);
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
