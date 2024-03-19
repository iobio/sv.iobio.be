/**
 * Function will take a url and return an array of chromosomes similar to the bands function. Chromosomes are not zipped 
 * because they are small enough.
 */

import fs from 'fs/promises';

async function parseChromosomes(url) {
    const fileBuffer = await fs.readFile(url);

    const data = fileBuffer.toString('utf-8');
    return processChromosomes(data);
}

function processChromosomes(dataString) {
    // Assumes dataString is a string containing the entire data
    const lines = dataString.split("\n");
    const chromosomes = lines.map(line => {
        const parts = line.split("\t");
        /**
         * These files have 9 columns, the last of which contains semi-colon, then equal sign, separated values with the chromosome name in there.
         * Index 4 is the length of the chromosome. Index 8 is the last column.
         */

        if (parts.length === 9) {
            const infoRow = parts[8].split(';');
            //turn into a dictionary where the key is the first part and the value is the second part after the equal sign
            const info = infoRow.reduce((acc, item) => {
                const [key, value] = item.split('=');
                acc[key] = value;
                return acc;
            }, {});

            const chromosomeName = info['Name'];

            return {
                chr: chromosomeName,
                length: parseInt(parts[4], 10),
            };
        }
        
    }).filter(chromosome => chromosome); // Remove any undefined entries due to filtering
    return chromosomes;
}

export { parseChromosomes };