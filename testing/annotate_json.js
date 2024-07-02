/**
 * This is for testing only
 * 
 * This works on data files we have in /Users/emerson/Documents/Data/SV.iobio_testData/svpipe_results/Smoove || /Users/emerson/Documents/Data/SV.iobio_testData/svpipe_results/Manta
 * 
 * 
 */

import fs from 'fs';
import path from 'path';
import { execSync, spawn } from 'child_process';

// const smooveDataLoc = '/Users/emerson/Documents/Data/SV.iobio_testData/svpipe_results/Smoove';
const mantaDataLoc = '/Users/emerson/Documents/Data/SV.iobio_testData/svpipe_results/Manta';

function annotateJson() {
    //get the json file
    const dataDir = mantaDataLoc;
    const files = fs.readdirSync(dataDir);
    const jsonFile = files.filter(file => file.endsWith('.json'));

    //there should only be one so we can just grab the first one
    const jsonPath = path.join(dataDir, jsonFile[0]);
    //get the name of the json file just take off the .json at the end
    const fileName = jsonFile[0].slice(0, -5);
    

    //read the json file
    const data = fs.readFileSync(jsonPath);

    //parse the json file
    const jsonData = JSON.parse(data);

    let outputJson = [];
    let vcfPath = path.join(dataDir, fileName + '.vcf.gz');
    let locationMap = {};

    //For each item in the json file we want to add the vcf annotation to the json file using bcftools and the vcf file
    for (const [index, sample] of jsonData.entries()) {
        let variantInfo = sample.variantEvaluations[0];
        let variantLocation = `chr${variantInfo.contigName}:${variantInfo.start}-${variantInfo.end}`;

        //get the vcf annotation for the variant
        let bcftoolsCmd = `bcftools view -r ${variantLocation} -H ${vcfPath} 2>/dev/null`;
        let vcfOutput = execSync(bcftoolsCmd).toString().split('\t');
        let vcfInfo = vcfOutput[7];
        let rank = index + 1;
        let gene = { geneSymbol: sample.geneSymbol, geneId: sample.geneIdentifier.geneId, ensemblId: sample.geneIdentifier.ensemblId};
        let exomiserCombScore = sample.combinedScore;
        let exomiserPval = sample.pValue;
        let exomiserPriorityScore = sample.priorityScore;

        //add the vcf annotation to the variant info 
        variantInfo.vcfInfo = vcfInfo;
        variantInfo.rank = rank;
        variantInfo.gene = gene;
        variantInfo.exomiserCombScore = exomiserCombScore;
        variantInfo.exomiserPval = exomiserPval;
        variantInfo.exomiserPriorityScore = exomiserPriorityScore;
        variantInfo.otherGenes = [];

        //if we havent seen this location before we want to add it to the location map
        if (!locationMap[variantLocation]) {
            locationMap[variantLocation] = variantInfo;
        } else {
            //if we have seen this location so we will add this variant info to the otherGenes array
            locationMap[variantLocation].otherGenes.push(gene);
        }
    }

    //add the location map to the output json
    for (const [location, variantInfo] of Object.entries(locationMap)) {
        outputJson.push(variantInfo);
    }

    //ultimately we will return this annotated json file
    return outputJson;
}

function vcfToJson(filePath, callback) {
    console.log(filePath);
    // Determine if the filePath is a URL
    const isUrl = filePath.startsWith('http://') || filePath.startsWith('https://') || filePath.startsWith('ftp://');

    let bcftoolsCmd;
    if (isUrl) {
        // Use curl to stream the file from the URL
        bcftoolsCmd = spawn('sh', ['-c', `curl -s -k ${filePath} | bcftools view -H -`]);
    } else {
        // Directly use bcftools for local files
        bcftoolsCmd = spawn('bcftools', ['view', '-H', filePath]);
    }

    let outputJson = [];
    let buffer = '';

    bcftoolsCmd.stdout.on('data', (data) => {
        buffer += data.toString();
        let lines = buffer.split('\n');
        buffer = lines.pop(); // Save the incomplete line
        
        let validChrom = new Set([
            ...Array.from({length:22}, (_, i) => (i + 1).toString()), 
            'X', 'Y'
        ]);

        for (const line of lines) {
            if (!line.trim()) continue;

            let variant = line.split('\t');
            let infoFields = variant[7].split(';');
            let end = variant[1];
            for (let field of infoFields) {
                if (field.startsWith('END=')) {
                    end = field.split('=')[1];
                    break;
                }
            }

            let type = 'none';
            for (let field of infoFields) {
                if (field.startsWith('SVTYPE=')) {
                    type = field.split('=')[1];
                    break;
                }
            }

            let variantInfo = {
                contigName: variant[0].replace(/^chr/, ''),
                start: parseInt(variant[1]),
                end: parseInt(end),
                ref: variant[3],
                alt: variant[4],
                //for now match
                variantLocation:`chr${variant[0].replace(/^chr/, '')}:${variant[1]}-${end}`, 
                vcfInfo: '',
                rank: '',
                gene: '',
                exomiserCombScore: 0,
                exomiserPval: 0,
                exomiserPriorityScore: 0,
                otherGenes: [], 
                type: type
            };

            if (!validChrom.has(variantInfo.contigName)){
                continue;
            }

            //filter make sure the variant is at least 100bp
            if (variantInfo.end - variantInfo.start <= 500) {
                continue;
            }

            outputJson.push(variantInfo);
        }
    });

    bcftoolsCmd.stderr.on('data', (data) => {
        console.error(`stderr: ${data}`);
    });

    bcftoolsCmd.on('close', (code) => {
        if (buffer.trim()) {
            let variant = buffer.split('\t');
            let infoFields = variant[7].split(';');
            let end = variant[1];
            for (let field of infoFields) {
                if (field.startsWith('END=')) {
                    end = field.split('=')[1];
                    break;
                }
            }

            let variantInfo = {
                contigName: variant[0],
                start: variant[1],
                end: end,
                ref: variant[3],
                alt: variant[4]
            };
            outputJson.push(variantInfo);
        }

        callback(outputJson);
    });
}

export { annotateJson, vcfToJson};
