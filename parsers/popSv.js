/**
 * Use bcftools to grab SVs from the SVAFotate_popAFs_GRCh38.sorted.bed.gz files for a given chromosome.
 */

import { spawn } from "child_process";

const popSvFile38 = "data/SVAFotate_popAFs_GRCh38.sorted.v4.1.bed.gz";
const popSvIndex38 = "data/SVAFotate_popAFs_GRCh38.sorted.v4.1.bed.gz.tbi";
const popSvFile37 = "data/SVAFotate_popAFs_GRCh37.sorted.v4.1.bed.gz";
const popSvIndex37 = "data/SVAFotate_popAFs_GRCh37.sorted.v4.1.bed.gz.tbi";

async function grabChrPopSvs(chr, start, end, svlen, build, prefix) {
    if (chr.startsWith("chr")) {
        chr = chr.slice(3);
    }
    if (build !== "hg38" && build !== "hg19") {
        throw new Error("Invalid build. Use 'GRCh38' or 'GRCh37'.");
    }

    const popSvFile = build === "hg38" ? popSvFile38 : popSvFile37;
    const popSvFileFullPath = `${prefix}${popSvFile}`;

    const s = Math.max(0, start - (parseInt(svlen) + 2)); // 2 is a small buffer to ensure we get all the SVs that are around
    const e = end + (parseInt(svlen) + 2); // 2 is a small buffer to ensure we get all the SVs that are around

    // use tabix to grab the SVs from the bed file
    const tabixCmd = spawn("tabix", [popSvFileFullPath, `${chr}:${s}-${e}`]);

    let buffer = "";
    let errorBuffer = "";
    let output = [];

    return new Promise((resolve, reject) => {
        tabixCmd.stdout.on("data", (data) => {
            buffer += data.toString();
            let lines = buffer.split("\n");
            buffer = lines.pop(); // Save the incomplete line

            let svs = formatSvs(lines);
            output = output.concat(svs);
        });

        tabixCmd.stderr.on("data", (data) => {
            errorBuffer += data.toString();
        });

        tabixCmd.on("close", (code) => {
            if (code === 0) {
                resolve(output);
            } else {
                reject(`Tabix failed with code ${code}`);
            }
        });
    });
}

function formatSvs(lines) {
    const svs = lines
        .map((line) => {
            const parts = line.split("\t");
            if (parts.length >= 173) {
                return {
                    chr: parts[0],
                    start: parseInt(parts[1], 10),
                    end: parseInt(parts[2], 10),
                    svlen: parts[3],
                    svtype: parts[4],
                    source: parts[5],
                    sv_id: parts[6],
                    af: parts[7],
                    num_hom_ref: parts[8],
                    num_het: parts[9],
                    num_hom_alt: parts[10],
                    pop_max_af: parts[172],
                };
            }
        })
        .filter((sv) => sv); // Remove any undefined entries due to filtering
    return svs;
}

export { grabChrPopSvs };
