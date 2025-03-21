import { parseBands } from "./parsers/bands.js";
import { parseChromosomes } from "./parsers/chromosomes.js";
import { parseClinGenRegions, parseClinGenGenes } from "./parsers/clinGen.js";
import { parseHg19Centromeres, parseHg38Centromeres } from "./parsers/centromeres.js";
import { vcfToJson, vcfSamples, vcfQuality, getPopSvs } from "./testing/annotate_json.js";
import { getOverlappedGenes, getGeneAssociations, getCanonicalTranscript } from "./testing/dbHelpers.js";
import sqlite3 from "sqlite3";
import express from "express";

const app = express();
const port = 7477;

// let prefix = "./"; //development
let prefix = "/"; //production

//will need to set the cors headers to allow all
app.use((req, res, next) => {
    res.header("Access-Control-Allow-Origin", "*");
    res.header("Access-Control-Allow-Headers", "Origin, X-Requested-With, Content-Type, Accept");
    next();
});

//json
app.use(express.json());

app.get("/", (req, res) => {
    res.send("sv.iobio Backend Running...");
});

//the bands endpoint expects a build hg19 or hg38 returns the bands
//The gzipped files are in our data folder
app.get("/bands", async (req, res) => {
    const build = req.query.build;
    let url;

    if (!build || (build !== "hg19" && build !== "hg38")) {
        res.status(400).send("valid build query parameter is required");
    } else {
        if (build === "hg19") {
            url = `${prefix}data/cytoBand_hg19.txt.gz`;
        } else {
            url = `${prefix}data/cytoBand_hg38.txt.gz`;
        }

        //Use the parseBands function to get the bands
        try {
            const bands = await parseBands(url);
            res.send(bands);
        } catch (e) {
            res.status(500).send(e.message);
        }
    }
});

//the chromosomes endpoint expects a build hg19 or hg38 returns the chromosomes
//The files are in our data folder and are not gzipped
app.get("/chromosomes", async (req, res) => {
    const build = req.query.build;
    let url;

    if (!build || (build !== "hg19" && build !== "hg38")) {
        res.status(400).send("valid build query parameter is required");
    } else {
        if (build === "hg19") {
            url = `${prefix}data/chromosomes_hg19.txt`;
        } else {
            url = `${prefix}data/chromosomes_hg38.txt`;
        }

        //Use the parseChromosomes function to get the chromosomes
        try {
            const chromosomes = await parseChromosomes(url);
            res.send(chromosomes);
        } catch (e) {
            res.status(500).send(e.message);
        }
    }
});

//the centromeres endpoint expects a build hg19 or hg38 returns the centromeres
//The files are in our data folder
app.get("/centromeres", async (req, res) => {
    const build = req.query.build;
    let url;

    if (!build || (build !== "hg19" && build !== "hg38")) {
        res.status(400).send("valid build query parameter is required");
    } else {
        if (build === "hg19") {
            url = `${prefix}data/gaps_ref_hg19.txt.gz`;
        } else {
            url = `${prefix}data/centromeres_hg38.txt`;
        }

        //Use the parseHg19Centromeres or parseHg38Centromeres function to get the centromeres
        try {
            let centromeres;
            if (build === "hg19") {
                centromeres = await parseHg19Centromeres(url);
            } else {
                centromeres = await parseHg38Centromeres(url);
            }
            res.send(centromeres);
        } catch (e) {
            res.status(500).send(e.message);
        }
    }
});

app.get("/genes", (req, res) => {
    let build = req.query.build;
    let source = req.query.source;
    let sourceText = "";

    if (!build || (build !== "hg19" && build !== "hg38")) {
        res.status(400).send("Valid build query parameter is required");
        return;
    }

    if (source === "refseq") {
        sourceText = ' AND source = "refseq"';
    }

    let query = "";
    if (build === "hg19") {
        query = 'SELECT * FROM genes WHERE build = "GRCh37"' + sourceText; //this is okay because I made this here it didnt come from the request
    } else if (build === "hg38") {
        query = 'SELECT * FROM genes WHERE build = "GRCh38"' + sourceText;
    }

    const db = new sqlite3.Database(`${prefix}data/geneinfo.db/gene.iobio.db`);

    db.all(query, [], (err, rows) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }
        //for each row we want to have this become a json object where the gene_symbol is the key
        let geneMap = {};
        rows.forEach((row) => {
            geneMap[row.gene_symbol] = row;
        });
        res.send(geneMap);
    });
});

app.get("/doseSensitiveGenes", async (req, res) => {
    let build = req.query.build;

    if (!build || (build !== "hg19" && build !== "hg38")) {
        res.status(400).send("Valid build query parameter is required");
        return;
    }

    const urlStart = `${prefix}/data/ClinGen_gene_`;

    try {
        let sensitiveGenes;
        sensitiveGenes = await parseClinGenGenes(build, urlStart);

        let response = { sensitiveGenes };

        res.send(response);
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.get("/doseSensitiveRegions", async (req, res) => {
    let build = req.query.build;

    if (!build || (build !== "hg19" && build !== "hg38")) {
        res.status(400).send("Valid build query parameter is required");
        return;
    }

    const urlStart = `${prefix}/data/ClinGen_region_`;

    try {
        let sensitiveRegions;
        sensitiveRegions = await parseClinGenRegions(build, urlStart);

        let response = { sensitiveRegions };

        res.send(response);
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.get("/genes/region", async (req, res) => {
    let build = req.query.build;
    let source = req.query.source;
    let startChr = req.query.startChr;
    let startPos = req.query.startPos;
    let endChr = req.query.endChr;
    let endPos = req.query.endPos;
    let sourceText = "";

    if (!build | !source | !startChr | !startPos | !endChr | !endPos) {
        res.status(400).send(
            "Endpoint requires a start chr & position as well as an end chr & position. Typical build and source are also required",
        );
        return;
    }

    const db = new sqlite3.Database(`${prefix}data/geneinfo.db/gene.iobio.db`);

    let geneMap = await getOverlappedGenes(build, source, startChr, startPos, endChr, endPos, sourceText, db);
    db.close();
    res.send(geneMap);
});

app.get("/sv/pop_overlaps", async (req, res) => {
    let build = req.query.build;
    let chr = req.query.chr;
    let start = req.query.start;
    let end = req.query.end;
    let svlen = req.query.svlen;

    if (!build || !chr || !start || !end) {
        res.status(400).send("Endpoint requires build, chr, start, and end query parameters");
        return;
    }

    if (build !== "hg19" && build !== "hg38") {
        res.status(400).send("Valid build query parameter is required");
        return;
    }

    let svs;
    try {
        svs = await getPopSvs(chr, start, end, svlen, build, prefix);
        res.send(svs);
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.post("/sv/info/batch", async (req, res) => {
    let build = req.query.build;
    let source = req.query.source;

    //the request body will have the variants: which is an array of objects with chromosome, startPos, endPos
    let variants = req.body.variants;

    if (!build | !source) {
        res.status(400).send("Endpoint requires typical build and source parameters");
        return;
    }

    const geneDb = new sqlite3.Database(`${prefix}data/geneinfo.db/gene.iobio.db`);
    const hpoDb = new sqlite3.Database(`${prefix}data/hpo.db`);

    let mappedVariants = [];

    for (let variant of variants) {
        let chromosome = `chr${variant.chromosome}`;
        let geneMap = await getOverlappedGenes(build, source, chromosome, variant.start, chromosome, variant.end, source, geneDb);
        let genes = Object.keys(geneMap);
        let { phenToGene, diseaseToGene } = await getGeneAssociations(genes, hpoDb);

        for (let gene_name of genes) {
            if (!geneMap.hasOwnProperty(gene_name)) {
                continue;
            }
            geneMap[gene_name].phenotypes = phenToGene[gene_name] || {};
            geneMap[gene_name].diseases = diseaseToGene[gene_name] || {};
        }
        variant.overlappedGenes = geneMap;
        mappedVariants.push(variant);
    }

    geneDb.close();
    hpoDb.close();
    res.send(mappedVariants);
});

//We will also want to get phenotypes for a given gene
app.get("/geneAssociations", async (req, res) => {
    let genes = req.query.genes ? req.query.genes.split(",") : [];

    if (!genes) {
        res.status(400).send("A list of genes is a required parameter for this endpoint");
        return;
    }

    const db = new sqlite3.Database(`${prefix}data/hpo.db`);
    const { phenToGene, diseaseToGene } = await getGeneAssociations(genes, db);
    db.close();

    res.send({ phenToGene, diseaseToGene });
});

app.get("/transcripts", (req, res) => {
    let genes = req.query.genes ? req.query.genes.split(",") : [];
    let build = req.query.build;
    let source = req.query.source;

    if (!genes.length || !build || !source) {
        res.status(400).send("A list of genes, a source, and a build are required parameters for this endpoint");
        return;
    }

    if (build == "hg19") {
        build = "GRCh37";
    } else if (build == "hg38") {
        build = "GRCh38";
    } else {
        res.status(400).send("Valid build query parameter is required");
        return;
    }

    const db = new sqlite3.Database(`${prefix}data/geneinfo.db/gene.iobio.db`);

    let placeholders = genes.map(() => "?").join(",");
    let query = `SELECT * FROM transcripts WHERE gene_name IN (${placeholders}) AND source = ? AND build = ?`;
    let params = [...genes, source, build];

    db.all(query, params, (err, rows) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }

        let transcriptMap = {};
        genes.forEach((gene) => {
            let geneRows = rows.filter((row) => row.gene_name === gene);
            let canonical = getCanonicalTranscript(geneRows);
            transcriptMap[gene] = canonical;
        });
        res.send(transcriptMap);
    });
});

app.get("/phenotypeGenes", (req, res) => {
    let phenotypes = req.query.phenotypes ? req.query.phenotypes.split(",") : [];

    if (!phenotypes) {
        res.status(400).send("A list of phenotypes is a required parameter for this endpoint");
        return;
    }

    const db = new sqlite3.Database(`${prefix}data/hpo.db`);

    let placeholders = phenotypes.map(() => "?").join(",");

    let geneQuery = `SELECT term_to_gene.*, genes.gene_symbol FROM term_to_gene JOIN genes ON term_to_gene.gene_id = genes.gene_id WHERE term_to_gene.term_id IN (${placeholders})`;

    db.all(geneQuery, phenotypes, (err, rows) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }

        let geneMap = {};
        rows.forEach((row) => {
            if (geneMap.hasOwnProperty(row.gene_symbol)) {
                geneMap[row.gene_symbol][row.term_id] = row;
            } else {
                geneMap[row.gene_symbol] = {};
                geneMap[row.gene_symbol][row.term_id] = row;
            }
        });
        res.send(geneMap);
    });
});

app.get("/phenotype/lookup", (req, res) => {
    let phenotype = req.query.phenotype;
    let hpo_id = req.query.hpo_id;

    if (!phenotype && !hpo_id) {
        res.status(400).send("A phenotype or hpo_id is a required parameter for this endpoint");
        return;
    } else if (phenotype && hpo_id) {
        res.status(400).send("Only one of phenotype or hpo_id is allowed");
        return;
    }

    let query;
    let params;
    const db = new sqlite3.Database(`${prefix}data/hpo.db`);

    if (phenotype) {
        query = `SELECT term_id, name FROM Terms WHERE name LIKE ?`;
        params = `%${phenotype}%`;
    } else {
        query = `SELECT term_id, name FROM Terms WHERE term_id = ?`;
        params = hpo_id;
    }

    db.all(query, params, (err, row) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }

        res.send(row);
    });
});

app.get("/dataFromVcf", async (req, res) => {
    let vcfPath = req.query.vcfPath;
    let sampleName = req.query.sampleName;
    let build = req.query.build;

    let bandUrl;

    if (build === "hg19") {
        bandUrl = `${prefix}data/cytoBand_hg19.txt.gz`;
    } else {
        bandUrl = `${prefix}data/cytoBand_hg38.txt.gz`;
    }

    let bandList = null;
    try {
        bandList = await parseBands(bandUrl);
    } catch (e) {
        console.error("Error parsing bands file: ", e);
    }

    if (!vcfPath) {
        res.status(400).send("Valid vcfPath query parameter is required");
        return;
    }

    let json;
    try {
        if (!sampleName) {
            json = vcfToJson(
                vcfPath,
                (jsonOutput) => {
                    res.send(jsonOutput);
                },
                bandList,
            );
        } else {
            json = vcfToJson(
                vcfPath,
                (jsonOutput) => {
                    res.send(jsonOutput);
                },
                bandList,
                sampleName,
            );
        }
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.get("/vcfSamples", async (req, res) => {
    let vcfPath = req.query.vcfPath;

    if (!vcfPath) {
        res.status(400).send("Valid vcfPath query parameter is required");
        return;
    }

    let json;
    try {
        json = vcfSamples(vcfPath, (jsonOutput) => {
            res.send(jsonOutput);
        });
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.get("/vcfQuality", async (req, res) => {
    let vcfPath = req.query.vcfPath;
    let sampleName = req.query.sampleName;

    if (!vcfPath) {
        res.status(400).send("Valid vcfPath query parameter is required");
        return;
    }

    let json;
    try {
        if (sampleName) {
            json = vcfQuality(
                vcfPath,
                (jsonOutput) => {
                    res.send(jsonOutput);
                },
                sampleName,
            );
        } else {
            json = vcfQuality(vcfPath, (jsonOutput) => {
                res.send(jsonOutput);
            });
        }
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.listen(port, () => {
    console.log(`Example app listening at http://localhost:${port}`);
});

export default app;
