import { parseBands } from "./parsers/bands.js";
import express from 'express';

const app = express();
const port = 3000;

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

app.listen(port, () => {
    console.log(`Example app listening at http://localhost:${port}`);
});

export default app;