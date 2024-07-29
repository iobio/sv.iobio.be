import { spawn } from 'child_process';

function vcfToJson(filePath, callback, sampleName=null) {
    // Determine if the filePath is a URL
    const isUrl = filePath.startsWith('http://') || filePath.startsWith('https://') || filePath.startsWith('ftp://');

    let bcftoolsCmd;
    if (isUrl) {
        if (!sampleName) {
            // Use curl to stream the file from the URL
            bcftoolsCmd = spawn('sh', ['-c', `curl -s -k ${filePath} | bcftools view -H`]);            
        } else {
            // Use curl to stream the file from the URL
            bcftoolsCmd = spawn('sh', ['-c', `curl -s -k ${filePath} | bcftools view --samples ${sampleName} -H`]);
        }
    } else {
        if (!sampleName) {
            // Directly use bcftools for local files
            bcftoolsCmd = spawn('bcftools', ['view', '-H', filePath]);  
        } else {
            // Directly use bcftools for local files
            bcftoolsCmd = spawn('bcftools', ['view', '-H', '--samples', sampleName, filePath]);
        }
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
            //if variant @ 9 starts with 0/0 skip it because it is a reference
            if (variant[9].startsWith('0/0')) {
                continue;
            }

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
                variantLocation:`chr${variant[0].replace(/^chr/, '')}:${variant[1]}-${end}`, 
                vcfInfo: '',
                overlappedGenes: {},
                type: type,
                genotype: variant[9]
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
        callback(outputJson);
    });
}

function vcfSamples(filePath, callback) {
    // Determine if the filePath is a URL
    const isUrl = filePath.startsWith('http://') || filePath.startsWith('https://') || filePath.startsWith('ftp://');

    let bcftoolsCmd;
    if (isUrl) {
        // Use curl to stream the file from the URL
        bcftoolsCmd = spawn('sh', ['-c', `curl -s -k ${filePath} | bcftools query --list-samples`]);
    } else {
        // Directly use bcftools for local files
        bcftoolsCmd = spawn('sh', ['-c', `bcftools query --list-samples ${filePath}`]);
    }

    let outputJson = [];
    let buffer = '';

    bcftoolsCmd.stdout.on('data', (data) => {
        buffer += data.toString();
        let lines = buffer.split('\n');
        buffer = lines.pop(); // Save the incomplete line
        
        for (const line of lines) {
            if (!line.trim()) continue;

            outputJson.push(line);
        }
    });

    bcftoolsCmd.stderr.on('data', (data) => {
        console.error(`stderr: ${data}`);
    });

    bcftoolsCmd.on('close', (code) => {
        callback(outputJson);
    });
}

export { vcfToJson, vcfSamples };
