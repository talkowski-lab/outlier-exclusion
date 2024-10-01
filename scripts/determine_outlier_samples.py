import argparse
import duckdb


def write_sv_filter_outliers(con, filter_id, iqr_mult):
    if iqr_mult < 0:
        raise ValueError('IQR multiplier must be greater than or equal to 0')

    sql = f'SELECT quantile_cont(count, [0.25, 0.5, 0.75]) AS iqr FROM sv_counts_{filter_id};'
    iqr = con.sql(sql).fetchall()[0][0]
    sql = (f'CREATE OR REPLACE TABLE outliers_{filter_id}'
            ' AS'
            ' SELECT sample'
           f' FROM sv_counts_{filter_id}'
            ' WHERE count < $1 - $2 * $4 OR count > $1 + $3 * $4')
    con.execute(sql, [iqr[1], iqr[0], iqr[2], iqr_mult])


def find_sv_count_outliers(db, iqr_mult):
    with duckdb.connect(db) as con:
        filter_ids = con.sql('SELECT id FROM sv_filters;').fetchall()
        for i in filter_ids:
            write_sv_filter_outliers(con, i[0], iqr_mult)


def find_wgd_outliers(db, wgd_scores_path, min_wgd, max_wgd):
    if min_wgd > max_wgd:
        raise ValueError('Min WGD score must be >= max wgd score')

    with duckdb.connect(db) as con:
        sql = 'CREATE OR REPLACE TABLE wgd_scores (sample VARCHAR, score FLOAT);'
        con.sql(sql)
        sql = (f'COPY wgd_scores FROM \'{wgd_scores_path}\''
               '(FORMAT CSV, HEADER false, DELIMITER \'\t\')')
        con.sql(sql)
        sql = ('CREATE OR REPLACE TABLE wgd_outliers'
               ' AS SELECT sample'
               ' FROM wgd_scores'
               ' WHERE score < ? OR score > ?')
        con.execute(sql, [min_wgd, max_wgd])


def main(args):
    counts_db = Path(args.sv_counts_db)
    if not counts_db.is_file():
        raise FileNotFoundError('Counts database not found')
    wgd_scores_path = Path(args.wgd_scores)
    if not wgd_scores_path.is_file():
        raise FileNotFoundError('WGD scores file not found')

    find_sv_count_outliers(counts_db, args.iqr_mult)
    find_wgd_outliers(counts_db, wgd_scores_path, args.min_wgd, args.max_wgd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determine sample outliers in GATK-SV callset"
    )
    parser.add_argument('sv_counts_db',
                        metavar = 'SV_COUNTS_DB',
                        help ='Path to the input SV counts DuckDB database')
    parser.add_argument('iqr_mult',
                        metavar = 'IQR_MULTIPLIER',
                        help = 'SVs per genome IQR multiplier',
                        type = float)
    parser.add_argument('wgd_scores',
                        metavar = 'WGD_SCORES',
                        help = 'Path to the sample WGD scores')
    parser.add_argument('min_wgd',
                        metavar = 'MIN_WGD',
                        help = 'Minimum WGD score')
    parser.add_argument('max_wgd',
                        metavar = 'MAX_WGD',
                        help = 'Maximum WGD score')
    args = parser.parse_args()
