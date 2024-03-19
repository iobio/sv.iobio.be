/**
 * Function will take a url return an array of centromeres.
 * 
 * hg19 centromeres are in a gzipped file in the data folder and are in a different format than hg38 centromeres.
 * 
 * There will be a function for each build due to the different formats.
 * 
 */

import fs from 'fs/promises';
import zlib from 'zlib';
import { promisify } from 'util';

const gunzip = promisify(zlib.gunzip);

async function parseHg19Centromeres(url) {
    const fileBuffer = await fs.readFile(url);
    
    const decompressedBuffer = await gunzip(fileBuffer);

    const data = decompressedBuffer.toString('utf-8');
    return processHg19Centromeres(data);
}

function processHg19Centromeres(dataString) {
    const lines = dataString.split("\n");

    const centromeres = lines.map(line => {
        const parts = line.split("\t");
        
        // If the row is not a centromere skip it
        if (parts[7] !== "centromere") {
            return;
        } else {
            return {
                chr: parts[1],
                start: parseInt(parts[2], 10),
                end: parseInt(parts[3], 10),
            };
        }
    }).filter(centromere => centromere); // Remove any undefined entries due to filtering
    return centromeres;
}

async function parseHg38Centromeres(url) {
    const fileBuffer = await fs.readFile(url);

    const data = fileBuffer.toString('utf-8');
    return processHg38Centromeres(data);
}

function processHg38Centromeres(dataString) {
    const lines = dataString.split("\n");

    //Set up some variables to keep track of the current chromosome and the start and end of the centromere region
    let currentChromosome = null;
    let start = 0;
    let end = 0;

    const centromeres = lines.map(line => {
        const parts = line.split("\t");
        //if the currentChromosome is null, then we are at the beginning of the file and we need to set the currentChromosome and the start of the centromere region
        if (currentChromosome === null) {
            currentChromosome = parts[1];
            start = parseInt(parts[2], 10);
            end = parseInt(parts[3], 10);
        } else if (currentChromosome === parts[1]) {
            end = parseInt(parts[3], 10);
        } else {
            const centromere = {
                chr: currentChromosome,
                start: start,
                end: end,
            };
            //reset the start and end of the centromere region
            start = parseInt(parts[2], 10);
            end = parseInt(parts[3], 10);
            currentChromosome = parts[1];
            return centromere;
        }
    }).filter(centromere => centromere); // Remove any undefined entries due to filtering
    return centromeres;
}

export { parseHg19Centromeres, parseHg38Centromeres };