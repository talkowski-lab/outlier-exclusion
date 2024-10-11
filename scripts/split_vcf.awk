# Split an indexed VCF into N VCFs
# Usage: awk -f split_vcf.awk VCF N [PREFIX]

BEGIN {
    vcf = ARGV[1]
    n_splits = ARGV[2]
    if (ARGC == 4) {
        prefix = ARGV[3]
    } else {
        prefix = "./"
    }

    prefix_dir = dirname(prefix)
    if (prefix_dir && system("test -d " prefix_dir) == 1) {
        if (system("mkdir -p " prefix_dir) != 0) {
            print_err("could not create directory " prefix_dir)
            exit 87
        }
    }

    cmd = "bcftools index --nrecords '" vcf "'"
    cmd | getline n_records
    close(cmd)

    n_records = int(n_records)
    if (n_records == 0) {
        print_err("VCF has 0 records")
        exit 86
    }

    if (n_splits > n_records) {
        n_splits = n_records
    }

    records_per_split = int(n_records / n_splits)
    max_counter = n_splits - 1
    counter_width = length(max_counter "")
    output_template = prefix "%0" counter_width "d_" basename(vcf)

    cmd = "bcftools head '" vcf "'"
    while ((cmd | getline line) > 0) {
        vcf_header[++i] = line
    }
    n_header_lines = i

    in_cmd = "bgzip --decompress --stdout '" vcf "'"
    i = -1
    j = records_per_split
    while ((in_cmd | getline line) > 0) {
        if (line ~ /^#/) {
            continue
        }

        if (i < n_splits - 1 && j++ == records_per_split) {
            if (out_cmd) {
                close(out_cmd)
            }

            output = sprintf(output_template, ++i)
            out_cmd = "bgzip --output '" output "'"
            j = 0

            for (k = 1; k <= n_header_lines; ++k) {
                print vcf_header[k] | out_cmd
            }
        }

        print line | out_cmd
    }

    if (out_cmd) {
        close(out_cmd)
    }
    close(in_cmd)
    exit 0
}

function basename(x) {
    sub(/^.*\//, "", x)

    return x
}

function dirname(x) {
    sub(/[^\/]+$/, "", x)

    return x
}

function print_err(x) {
    print x > "/dev/stderr"
}
