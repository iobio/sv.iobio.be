import { parseBands } from "./parsers/bands.js";
import { parseChromosomes } from "./parsers/chromosomes.js";
import { parseHg19Centromeres, parseHg38Centromeres } from "./parsers/centromeres.js";
import { annotateJson } from "./testing/annotate_json.js";
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

//the annotate endpoint is for testing only
app.get('/vcfjson', async (req, res) => {
    let annotatedJson;
    try {
        annotatedJson = await annotateJson();
        res.send(annotatedJson);
    } catch (e) {
        res.status(500).send(e.message);
    }
});

app.listen(port, () => {
    console.log(`Example app listening at http://localhost:${port}`);
});

export default app;