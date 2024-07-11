
//Get a list of ids given a list of names
export function getGeneIdsList(geneNamesList, db) {
    return new Promise((resolve, reject) => {
        let placeholders = geneNamesList.map(() => '?').join(',');
        let query = `SELECT gene_id FROM genes WHERE gene_symbol IN (${placeholders})`;

        db.all(query, geneNamesList, (err, rows) => {
            if (err) {
                reject(err);
            } else {
                let geneIdsList = rows.map(row => row.gene_id);
                resolve(geneIdsList);
            }
        })
    })
}

export function getGeneAssociations(genes, db) {
    /**
     * 
     */
    return new Promise((resolve, reject) => {
        getGeneIdsList(genes, db).then(geneIdsList => {
            let placeholders = geneIdsList.map(() => '?').join(',');

            const getPhenotypes = () => {
                return new Promise((resolve, reject) => {
                    let phenQuery = `SELECT term_to_gene.*, genes.gene_symbol, terms.name 
                        FROM term_to_gene 
                        JOIN genes ON term_to_gene.gene_id = genes.gene_id
                        JOIN terms ON term_to_gene.term_id = terms.term_id
                        WHERE term_to_gene.gene_id IN (${placeholders})`;

                    db.all(phenQuery, geneIdsList, (err, rows) => {
                        if (err) {
                            reject(err);
                        }
                
                        let phenotypesToGene = {};
                        //if thre are no rows then we want to return an empty object
                        if (!rows || rows.length === 0) {
                            //resolve an empty object
                        } else {
                            rows.forEach(row => {
                                if (phenotypesToGene.hasOwnProperty(row.gene_symbol)){
                                    phenotypesToGene[row.gene_symbol][row.term_id] = row;
                                } else {
                                    //establish the structure of the phenotypes to gene
                                    phenotypesToGene[row.gene_symbol] = {};
                                    phenotypesToGene[row.gene_symbol][row.term_id] = row;
                                }
                            })
                        }
                        resolve(phenotypesToGene);
                    });
                })
            };

            const getDiseases = () => {
                return new Promise((resolve, reject) => {
                    let diseaseQuery = `
                        SELECT gene_to_disease.*, genes.gene_symbol, diseases.disease_name 
                        FROM gene_to_disease 
                        JOIN genes ON gene_to_disease.gene_id = genes.gene_id
                        JOIN diseases ON gene_to_disease.disease_id = diseases.disease_id 
                        WHERE gene_to_disease.gene_id IN (${placeholders})`
    
                    db.all(diseaseQuery, geneIdsList, (err, rows) => {
                        if (err) {
                            reject(err);
                        }
                        let diseasesToGene = {};

                        //if there are no rows then we want to return an empty object
                        if (!rows || rows.length === 0) {
                            //resolve an empty object
                        } else {
                            rows.forEach(row => {
                                if (diseasesToGene.hasOwnProperty(row.gene_symbol)){
                                    diseasesToGene[row.gene_symbol][row.disease_id] = row;
                                } else {
                                    diseasesToGene[row.gene_symbol] = {};
                                    diseasesToGene[row.gene_symbol][row.disease_id] = row;
                                }
                            })
                        }
                        resolve(diseasesToGene);
                    });
                })
            };

            Promise.all([getPhenotypes(), getDiseases()]).then(([phenotypeMap, diseaseMap]) => {
                resolve({phenToGene: phenotypeMap, diseaseToGene: diseaseMap});
            }).catch(err => {
                reject(err);
            });
        });
    });
}

export function getOverlappedGenes(build, source, startChr, startPos, endChr, endPos, sourceText, db) {
    /**
     * 
     */
    return new Promise((resolve, reject) => {
        if (source === 'refseq') {
            sourceText = ' AND source = "refseq"';
        }
    
        let query = '';
        let buildText = '';
    
        if (build === 'hg19') {
            buildText = 'build = "GRCh37"';
        } else if (build === 'hg38') {
            buildText = 'build = "GRCh38"';
        }
    
        let params = [];
    
        if (startChr == endChr) {
            params = [startChr, startPos, endPos, startPos, startPos, endPos, endPos]

            query = `SELECT * FROM genes WHERE ${buildText} AND chr = ? AND ((start >= ? AND end <= ?) OR ((? >= start AND ? < end) OR (? >= start AND ? < end)))` + sourceText;
        } else {
            params = [startChr, startPos, startPos, endChr, endPos, endPos]
            query = `SELECT * FROM genes WHERE ${buildText} AND ((chr = ? AND (start >= ? OR end >= ?)) OR (chr = ? AND (end <= ? OR start <= end ?))` + sourceText;
        }
    
        db.all(query, params, (err, rows) => {
            if (err) {
                reject(err);
            }
            //for each row we want to have this become a json object where the gene_symbol is the key
            let geneMap = {};
            rows.forEach(row => {
                //if the gene_symbols is null then we dont want to include it
                if (row.gene_symbol === null) {
                    return;
                }
                geneMap[row.gene_symbol] = row;
            });
            resolve(geneMap);
        });
    });
}