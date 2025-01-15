import { spawn } from "child_process";
import { grabChrPopSvs } from "../parsers/popSv.js";

function vcfToJson(filePath, callback, bandList, sampleName = null) {
    // Determine if the filePath is a URL
    const isUrl = filePath.startsWith("http://") || filePath.startsWith("https://") || filePath.startsWith("ftp://");

    let bcftoolsCmd;
    if (isUrl) {
        if (!sampleName) {
            // Use curl to stream the file from the URL
            bcftoolsCmd = spawn("sh", ["-c", `curl -s -k ${filePath} | bcftools view -H | grep 'SVTYPE='`]);
        } else {
            // Use curl to stream the file from the URL
            bcftoolsCmd = spawn("sh", [
                "-c",
                `curl -s -k ${filePath} | bcftools view --samples ${sampleName} -H | grep 'SVTYPE='`,
            ]);
        }
    } else {
        if (!sampleName) {
            // Directly use bcftools for local files and add grep
            bcftoolsCmd = spawn("sh", ["-c", `bcftools view -H ${filePath} | grep 'SVTYPE='`]);
        } else {
            // Directly use bcftools for local files and add grep
            bcftoolsCmd = spawn("sh", ["-c", `bcftools view --samples ${sampleName} -H ${filePath} | grep 'SVTYPE='`]);
        }
    }

    let outputJson = [];
    let buffer = "";

    let uniqueVariants = {};

    bcftoolsCmd.stdout.on("data", (data) => {
        buffer += data.toString();
        let lines = buffer.split("\n");
        buffer = lines.pop(); // Save the incomplete line

        let validChrom = new Set([...Array.from({ length: 22 }, (_, i) => (i + 1).toString()), "X", "Y"]);

        for (const line of lines) {
            if (!line.trim()) continue;

            let variant = line.split("\t");
            //if variant @ 9 starts with 0/0 or ./. skip it because it is a reference
            if (variant[9].startsWith("0/0") || variant[9].startsWith("./.")) {
                continue;
            }

            let infoFields = variant[7].split(";");
            let end = variant[1];
            for (let field of infoFields) {
                if (field.startsWith("END=")) {
                    end = field.split("=")[1];
                    break;
                }
            }

            let type = false;
            let size = false;
            for (let field of infoFields) {
                if (field.startsWith("SVTYPE=")) {
                    type = field.split("=")[1];
                }
                if (field.startsWith("SVLEN=")) {
                    size = field.split("=")[1];
                }

                if (type && size) {
                    break;
                }
            }

            let variantInfo = {
                contigName: variant[0].replace(/^chr/, ""),
                start: parseInt(variant[1]),
                end: parseInt(end),
                size: size,
                quality: variant[5],
                variantLocation: `chr${variant[0].replace(/^chr/, "")}:${variant[1]}-${end}`,
                vcfInfo: "",
                overlappedGenes: {},
                type: type,
                genotype: variant[9],
                bands: [],
            };

            let key = `${variantInfo.contigName}:${variantInfo.start}-${variantInfo.end}-${variantInfo.type}`;
            if (uniqueVariants[key]) {
                continue;
            } else {
                uniqueVariants[key] = true;
            }

            if (!validChrom.has(variantInfo.contigName)) {
                continue;
            }

            if (bandList) {
                for (let band of bandList) {
                    if (band.chr === variant[0]) {
                        if (band.start <= variantInfo.start && band.end >= variantInfo.end) {
                            variantInfo.bands.push(band);
                        } else if (band.start >= variantInfo.start && band.start <= variantInfo.end) {
                            variantInfo.bands.push(band);
                        } else if (band.end >= variantInfo.start && band.end <= variantInfo.end) {
                            variantInfo.bands.push(band);
                        }
                    }
                }
            }

            //filter make sure the variant is at least 50bp
            if (size < 50 && variantInfo.end - variantInfo.start < 50) {
                continue;
            }

            outputJson.push(variantInfo);
        }
    });

    bcftoolsCmd.stderr.on("data", (data) => {
        console.error(`stderr: ${data}`);
    });

    bcftoolsCmd.on("close", (code) => {
        callback(outputJson);
    });
}

function vcfSamples(filePath, callback) {
    // Determine if the filePath is a URL
    const isUrl = filePath.startsWith("http://") || filePath.startsWith("https://") || filePath.startsWith("ftp://");

    let bcftoolsCmd;
    if (isUrl) {
        // Use curl to stream the file from the URL
        bcftoolsCmd = spawn("sh", ["-c", `curl -s -k ${filePath} | bcftools query --list-samples`]);
    } else {
        // Directly use bcftools for local files
        bcftoolsCmd = spawn("sh", ["-c", `bcftools query --list-samples ${filePath}`]);
    }

    let outputJson = [];
    let buffer = "";

    bcftoolsCmd.stdout.on("data", (data) => {
        buffer += data.toString();
        let lines = buffer.split("\n");
        buffer = lines.pop(); // Save the incomplete line

        for (const line of lines) {
            if (!line.trim()) continue;

            outputJson.push(line);
        }
    });

    bcftoolsCmd.stderr.on("data", (data) => {
        console.error(`stderr: ${data}`);
    });

    bcftoolsCmd.on("close", (code) => {
        callback(outputJson);
    });
}

function vcfQuality(filePath, callback, sampleName = null) {
    // Determine if the filePath is a URL
    const isUrl = filePath.startsWith("http://") || filePath.startsWith("https://") || filePath.startsWith("ftp://");

    let bcftoolsCmd;
    if (isUrl) {
        if (!sampleName) {
            // Use curl to stream the file from the URL
            bcftoolsCmd = spawn("sh", ["-c", `curl -s -k ${filePath} | bcftools query -f '%QUAL\n'`]);
        } else {
            // Use curl to stream the file from the URL
            bcftoolsCmd = spawn("sh", ["-c", `curl -s -k ${filePath} | bcftools query --samples ${sampleName} -f '%QUAL\n'`]);
        }
    } else {
        if (!sampleName) {
            // Directly use bcftools for local files
            bcftoolsCmd = spawn("sh", ["-c", `bcftools query -f '%QUAL\n' ${filePath}`]);
        } else {
            // Directly use bcftools for local files
            bcftoolsCmd = spawn("sh", ["-c", `bcftools query --samples ${sampleName} -f '%QUAL\n' ${filePath}`]);
        }
    }

    let qualityArray = [];
    let buffer = "";

    bcftoolsCmd.stdout.on("data", (data) => {
        buffer += data.toString();
        let lines = buffer.split("\n");
        buffer = lines.pop(); // Save the incomplete line

        for (const line of lines) {
            if (!line.trim()) continue;

            qualityArray.push(line);
        }
    });

    //Calculate the average quality, the median quality, and the IQR
    bcftoolsCmd.stderr.on("data", (data) => {
        console.error(`stderr: ${data}`);
    });

    let stats = {
        max: null,
        min: null,
        avg: null,
        median: null,
        stdDev: null,
    };

    bcftoolsCmd.on("close", (code) => {
        let quality = qualityArray.filter((value) => !isNaN(value) && value !== null && value !== undefined).map(Number);

        if (quality.length === 0) {
            callback(stats);
            return;
        }

        let sum = quality.reduce((a, b) => a + b, 0);
        let avg = sum / quality.length;
        let max = Math.max(...quality);
        let min = Math.min(...quality);

        quality.sort((a, b) => a - b);
        let median;
        if (quality.length % 2 == 0) {
            median = (quality[quality.length / 2 - 1] + quality[quality.length / 2]) / 2;
        } else {
            median = quality[(quality.length - 1) / 2];
        }

        //Calculate the standard deviation
        let sqDiffs = quality.map((value) => (value - avg) ** 2);
        let avgSqDiff = sqDiffs.reduce((a, b) => a + b, 0) / sqDiffs.length;
        let stdDev = Math.sqrt(avgSqDiff);

        stats = {
            max: max,
            min: min,
            avg: avg,
            median: median,
            stdDev: stdDev,
        };

        callback(stats);
    });
}

async function getPopSvs(chr, start, end, svlen, build, prefix) {
    /**
     * Returns all the population SVs that overlap with the given region
     *
     * Expects relative coordinates (coordinates based on the chromosome)
     */

    let chrSvs;
    const overlapProp = 0.8;
    let svs = [];

    try {
        chrSvs = await grabChrPopSvs(chr, start, end, svlen, build, prefix);

        if (!chrSvs) {
            console.log("No SVs found in the region");
            return [];
        }

        for (let sv of chrSvs) {
            if (
                (sv.start >= start && sv.start <= end) ||
                (sv.end >= start && sv.end <= end) ||
                (start >= sv.start && start <= sv.end) ||
                (end >= sv.start && end <= sv.end)
            ) {
                //one of the conditions for overlap is met
                let overlap = Math.min(end, sv.end) - Math.max(start, sv.start);
                let overlapProp1 = overlap / (end - start);
                let overlapProp2 = overlap / (sv.end - sv.start);

                //Reciprocal overlap similar to svafotate so if the other sv is giant or something like that it will not be included
                if (overlapProp1 >= overlapProp && overlapProp2 >= overlapProp) {
                    //add overlap proportion to sv
                    sv.overlapFractionProd = overlapProp1 * overlapProp2;
                    svs.push(sv);
                }
            }
        }
        return svs;
    } catch (error) {
        console.error(error);
        return [];
    }
}

export { vcfToJson, vcfSamples, vcfQuality, getPopSvs };
