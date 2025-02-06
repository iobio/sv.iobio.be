/**
 * Will parse the appropriate clingen file for now (grch37/38) genes and regions
 */
import fs from "fs/promises";

async function parseClinGenGenes(build, baseGeneUrl) {
    const fileUrl = baseGeneUrl + `${build}.tsv`;
    const fileBuffer = await fs.readFile(fileUrl);
    const fileString = fileBuffer.toString("utf-8");

    return _processClinGenGenes(fileString);
}

function _processClinGenGenes(fileString) {
    const lines = fileString.split("\n").slice(6);

    let sensitiveGenes = lines.reduce((acc, line, index) => {
        let ln = line.split("\t");
        //Key will be gene symbol because of how we use this downstream
        let k = ln[0];

        if (k == "") {
            return acc;
        }

        acc[k] = {
            geneSymbol: ln[0],
            geneId: ln[1],
            band: ln[2],
            location: ln[3],
            haploinsufficiency: {
                score: ln[4],
                description: ln[5],
                evidence: [ln[6], ln[7], ln[8], ln[9], ln[10], ln[11]],
                diseaseAssociation: ln[21],
            },
            triplosensitivity: {
                score: ln[12],
                description: ln[13],
                evidence: [ln[14], ln[15], ln[16], ln[17], ln[18], ln[19]],
                diseaseAssociation: ln[22],
            },
            updated: ln[20],
        };

        return acc;
    }, {});

    return sensitiveGenes;
}

async function parseClinGenRegions(build, baseRegionUrl) {
    const fileUrl = baseRegionUrl + `${build}.tsv`;
    const fileBuffer = await fs.readFile(fileUrl);
    const fileString = fileBuffer.toString("utf-8");

    return _processClinGenRegions(fileString);
}

function _processClinGenRegions(fileString) {
    const lines = fileString.split("\n").slice(6);

    let sensitiveRegions = lines.reduce((acc, line, index) => {
        let ln = line.split("\t");
        //Key will be the ISCA ID
        let k = ln[0];
        if (k == "") {
            return acc;
        }

        acc[k] = {
            iscaId: ln[0],
            iscaName: ln[1],
            band: ln[2],
            location: ln[3],
            haploinsufficiency: {
                score: ln[4],
                description: ln[5],
                evidence: [ln[6], ln[7], ln[8], ln[9], ln[10], ln[11]],
                diseaseAssociation: ln[21],
            },
            triplosensitivity: {
                score: ln[12],
                description: ln[13],
                evidence: [ln[14], ln[15], ln[16], ln[17], ln[18], ln[19]],
                diseaseAssociation: ln[22],
            },
            updated: ln[20],
        };

        return acc;
    }, {});

    return sensitiveRegions;
}

export { parseClinGenGenes, parseClinGenRegions };
