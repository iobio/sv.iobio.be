/**
 * For now I think one function exported that takes a URL (for a file) and returns an array of bands is sufficient.
 * Not making it "default" because I think we might want to add more functions to this file at a later time.
 * 
 * Bands have the following properties:
 * - chr (string)
 * - start (int)
 * - end (int) 
 * - name (string)
 * - gieStain (string)
 * 
 */

import fs from 'fs/promises';
import zlib from 'zlib';
import { promisify } from 'util';

const gunzip = promisify(zlib.gunzip);

async function parseBands(url) {
        const fileBuffer = await fs.readFile(url);

        const decompressedBuffer = await gunzip(fileBuffer);
        const data = decompressedBuffer.toString('utf-8');
        return processBands(data);
}

function processBands(dataString) {
    // Assumes dataString is a string containing the entire data
    const lines = dataString.split("\n");
    const bands = lines.map(line => {
        const parts = line.split("\t");
        if (parts.length === 5) {
            return {
                chr: parts[0],
                start: parseInt(parts[1], 10),
                end: parseInt(parts[2], 10),
                name: parts[3],
                gieStain: parts[4],
            };
        }
    }).filter(band => band); // Remove any undefined entries due to filtering
    return bands;
}

export { parseBands };