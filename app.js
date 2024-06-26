import { parseBands } from "./parsers/bands.js";
import { parseChromosomes } from "./parsers/chromosomes.js";
import { parseHg19Centromeres, parseHg38Centromeres } from "./parsers/centromeres.js";
import { annotateJson, vcfToJson } from "./testing/annotate_json.js";
import { getGeneIdsList } from "./testing/dbHelpers.js";
import sqlite3 from 'sqlite3';
import express from 'express';

const app = express();
const port = 3000;

//will need to set the cors headers to allow all
app.use((req, res, next) => {
    res.header('Access-Control-Allow-Origin', '*');
    res.header('Access-Control-Allow-Headers', 'Origin, X-Requested-With, Content-Type, Accept');
    next();
});

//the bands endpoint expects a build hg19 or hg38 returns the bands
//The gzipped files are in our data folder
app.get('/bands', async (req, res) => {
    const build = req.query.build;
    let url;

    if (!build || (build !== 'hg19' && build !== 'hg38')) {
        res.status(400).send('valid build query parameter is required');
    } else {
        if (build === 'hg19') {
            url = './data/cytoBand_hg19.txt.gz';
        } else {
            url = './data/cytoBand_hg38.txt.gz';
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
app.get('/chromosomes', async (req, res) => {
    const build = req.query.build;
    let url;

    if (!build || (build !== 'hg19' && build !== 'hg38')) {
        res.status(400).send('valid build query parameter is required');
    } else {
        if (build === 'hg19') {
            url = './data/chromosomes_hg19.txt';
        } else {
            url = './data/chromosomes_hg38.txt';
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
app.get('/centromeres', async (req, res) => {
    const build = req.query.build;
    let url;

    if (!build || (build !== 'hg19' && build !== 'hg38')) {
        res.status(400).send('valid build query parameter is required');
    } else {
        if (build === 'hg19') {
            url = './data/gaps_ref_hg19.txt.gz';
        } else {
            url = './data/centromeres_hg38.txt';
        }

        //Use the parseHg19Centromeres or parseHg38Centromeres function to get the centromeres
        try {
            let centromeres;
            if (build === 'hg19') {
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

app.get('/genes', (req, res) => {
    let build = req.query.build;
    let source = req.query.source;
    let sourceText = '';

    if (!build || (build !== 'hg19' && build !== 'hg38')) {
        res.status(400).send('Valid build query parameter is required');
        return;
    }

    if (source === 'refseq') {
        sourceText = ' AND source = "refseq"';
    }

    let query = '';
    if (build === 'hg19') {
        query = 'SELECT * FROM genes WHERE build = "GRCh37"' + sourceText; //this is okay because I made this here it didnt come from the request
    } else if (build === 'hg38') {
        query = 'SELECT * FROM genes WHERE build = "GRCh38"' + sourceText;
    }

    const db = new sqlite3.Database('./data/geneinfo.db/gene.iobio.db');

    db.all(query, [], (err, rows) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }
        //for each row we want to have this become a json object where the gene_symbol is the key
        let geneMap = {};
        rows.forEach(row => {
            geneMap[row.gene_symbol] = row;
        });
        res.send(geneMap);
    });
});

app.get('/genes/region', (req, res) => {
    let build = req.query.build;
    let source = req.query.source;
    let startChr = req.query.startChr
    let startPos = req.query.startPos
    let endChr = req.query.endChr
    let endPos = req.query.endPos
    let between = req.query.between ? req.query.between.split(','): [];
    let sourceText = ''

    if (!build | !source | !startChr | !startPos | !endChr | !endPos) {
        res.status(400).send('Endpoint requires a start chr & position as well as an end chr & position. Typical build and source are also required');
        return;
    }

    if (source === 'refseq') {
        sourceText = ' AND source = "refseq"';
    }

    let query = '';
    let buildText = '';

    if (build === 'hg19') {
        buildText = 'build = "GRCh37"';
    } else if (build === 'hg38') {
        buildText = 'build = "GRCh38"';
    }

    let params = [];

    if (startChr == endChr) {
        params = [startChr, startPos, endPos, startPos, startPos, endPos, endPos]
        query = `SELECT * FROM genes WHERE ${buildText} AND chr = ? AND ((start >= ? AND end <= ?) OR ((? >= start AND ? < end) OR (? >= start AND ? < end)))` + sourceText;
    } else if (between.length > 0) {
        let placeholders = between.map(() => '?').join(', ')
        params = [startChr, startPos, startPos, endChr, endPos, endPos, ...between]
        query = `SELECT * FROM genes WHERE ${buildText} AND ((chr = ? AND (start >= ? OR end >= ?)) OR (chr = ? AND (end <= ? OR start <= end ?)) OR (chr IN (${placeholders})))` + sourceText;   
    } else {
        params = [startChr, startPos, startPos, endChr, endPos, endPos]
        query = `SELECT * FROM genes WHERE ${buildText} AND ((chr = ? AND (start >= ? OR end >= ?)) OR (chr = ? AND (end <= ? OR start <= end ?))` + sourceText;
    }

    const db = new sqlite3.Database('./data/geneinfo.db/gene.iobio.db');

    db.all(query, params, (err, rows) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }
        //for each row we want to have this become a json object where the gene_symbol is the key
        let geneMap = {};
        rows.forEach(row => {
            geneMap[row.gene_symbol] = row;
        });
        res.send(geneMap);
    });
})

//We will also want to get phenotypes for a given gene
app.get('/genePhenotypes', (req, res) => {
    let geneName = req.query.gene;
    let geneNum = null;

    if (!geneName) {
        res.status(400).send('gene is a required parameter for this endpoint');
        return;
    }

    //We need to get the gene # from our genes table to use in other queries
    const db = new sqlite3.Database('./data/hpo.db');
    let geneQuery = 'SELECT * FROM genes WHERE gene_symbol = ?';
    let geneParams = [geneName];

    let phenQuery = ''

    db.all(geneQuery, geneParams, (err, rows) => {
        if (err) {
            res.status(500).send(err.message);
            return;
        }

        //There should only be one gene that matches a particular gene symbol if there aren't any dont continue
        if (rows.length == 0) {
            res.send([])
            return;
        }

        //If we get here we had a gene found with that symbol
        geneNum = rows[0].gene_id;
        phenQuery = `SELECT * FROM term_to_gene WHERE gene_id = ${geneNum}` //not washing because this came from our own db

        db.all(phenQuery, [], (err, rows) => {
            db.close();
    
            if (err) {
                res.status(500).send(err.message);
                return;
            }
    
            let phenMap = {}
            rows.forEach(row => {
                phenMap[row.term_id] = row;
            })
            res.send(phenMap)
        });
    })

})

//We will also want to get phenotypes for a given gene
app.get('/geneDiseases', (req, res) => {
    let geneName = req.query.gene;
    let geneNum = null;

    if (!geneName) {
        res.status(400).send('gene is a required parameter for this endpoint');
        return;
    }

    //We need to get the gene # from our genes table to use in other queries
    const db = new sqlite3.Database('./data/hpo.db');
    let geneQuery = 'SELECT * FROM genes WHERE gene_symbol = ?';
    let geneParams = [geneName];

    let diseaseQuery = ''

    db.all(geneQuery, geneParams, (err, rows) => {
        if (err) {
            res.status(500).send(err.message);
            return;
        }

        //There should only be one gene that matches a particular gene symbol if there aren't any dont continue
        if (rows.length == 0) {
            res.send([])
            return;
        }

        //If we get here we had a gene found with that symbol
        geneNum = rows[0].gene_id;
        diseaseQuery = `SELECT * FROM gene_to_disease WHERE gene_id = ${geneNum}`//not washing because this came from our own db

        db.all(diseaseQuery, [], (err, rows) => {
            db.close();
    
            if (err) {
                res.status(500).send(err.message);
                return;
            }
    
            let diseaseMap = {}
            rows.forEach(row => {
                diseaseMap[row.disease_id] = row;
            });
            res.send(diseaseMap)
        });
    })
})

//We will also want to get phenotypes for a given gene
app.get('/geneAssociations', async (req, res) => {
    let genes = req.query.genes ? req.query.genes.split(',') : [];

    if (!genes) {
        res.status(400).send('A list of genes is a required parameter for this endpoint');
        return;
    }

    const db = new sqlite3.Database('./data/hpo.db');
    let geneIdsList = await getGeneIdsList(genes, db);

    let placeholders = geneIdsList.map(() => '?').join(',');

    let phenQuery = `SELECT term_to_gene.*, genes.gene_symbol FROM term_to_gene JOIN genes ON term_to_gene.gene_id = genes.gene_id WHERE term_to_gene.gene_id IN (${placeholders})`;
    let phenPromise = new Promise((resolve, reject) => {
        db.all(phenQuery, geneIdsList, (err, rows) => {
            if (err) {
                reject(err);
            }
    
            let phenotypesToGene = {};
            rows.forEach(row => {
                if (phenotypesToGene.hasOwnProperty(row.gene_symbol)){
                    phenotypesToGene[row.gene_symbol][row.term_id] = row;
                } else {
                    //establish the structure of the phenotypes to gene
                    phenotypesToGene[row.gene_symbol] = {};
                    phenotypesToGene[row.gene_symbol][row.term_id] = row;
                }
            })
            resolve(phenotypesToGene);
        });
    })

    let diseaseQuery = `SELECT gene_to_disease.*, genes.gene_symbol FROM gene_to_disease JOIN genes ON gene_to_disease.gene_id = genes.gene_id WHERE gene_to_disease.gene_id IN (${placeholders})`
    let diseasePromise = new Promise((resolve, reject) => {
        db.all(diseaseQuery, geneIdsList, (err, rows) => {
            if (err) {
                reject(err);
            }
            let diseasesToGene = {};
            rows.forEach(row => {
                if (diseasesToGene.hasOwnProperty(row.gene_symbol)){
                    diseasesToGene[row.gene_symbol][row.disease_id] = row;
                } else {
                    diseasesToGene[row.gene_symbol] = {};
                    diseasesToGene[row.gene_symbol][row.disease_id] = row;
                }
            })
            resolve(diseasesToGene);
        });
    })
    let [phenotypeMap, diseaseMap] = await Promise.all([phenPromise, diseasePromise]);
    db.close();
    res.send({phenToGene: phenotypeMap, diseaseToGene: diseaseMap})
})

app.get('/phenotypeGenes', (req, res) => {
    let phenotypes = req.query.phenotypes ? req.query.phenotypes.split(',') : [];

    if (!phenotypes) {
        res.status(400).send('A list of phenotypes is a required parameter for this endpoint');
        return;
    }

    const db = new sqlite3.Database('./data/hpo.db');

    let placeholders = phenotypes.map(() => '?').join(',');

    let geneQuery = `SELECT term_to_gene.*, genes.gene_symbol FROM term_to_gene JOIN genes ON term_to_gene.gene_id = genes.gene_id WHERE term_to_gene.term_id IN (${placeholders})`;

    db.all(geneQuery, phenotypes, (err, rows) => {
        db.close(); // Close the database after fetching the data

        if (err) {
            res.status(500).send(err.message);
            return;
        }

        let geneMap = {};
        rows.forEach(row => {
            if (geneMap.hasOwnProperty(row.gene_symbol)){
                geneMap[row.gene_symbol][row.term_id] = row;
            } else {
                geneMap[row.gene_symbol] = {};
                geneMap[row.gene_symbol][row.term_id] = row;
            }
        });
        res.send(geneMap);
    });
});

//the annotate endpoint is for testing only
app.get('/vcfjson', async (req, res) => {
    let annotatedJson;
    try {
        annotatedJson = annotateJson();
        res.send(annotatedJson);
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.get('/dataFromVcf', (req, res) => {
    let vcfPath = req.query.vcfPath;

    if (!vcfPath) {
        res.status(400).send('Valid vcfPath query parameter is required');
        return;
    }

    let json;
    try {
        json = vcfToJson(vcfPath, (jsonOutput) => {
            res.send(jsonOutput);
        });
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.listen(port, () => {
    console.log(`Example app listening at http://localhost:${port}`);
});

export default app;