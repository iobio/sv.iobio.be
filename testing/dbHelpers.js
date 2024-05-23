import sqlite3 from 'sqlite3';

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
