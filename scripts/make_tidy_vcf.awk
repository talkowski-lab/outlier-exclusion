# Convert a VCF into tab-delimited format
# Usage: bgzip -cd VCF | awk -f make_tidy_vcf.awk [FILTERS] -
#
# FILTERS is a file of SV filters, one per line. Each filter is of the form:
# SVTYPE;MIN_SVLEN;MAX_SVLEN
# Only sites passing at least one of the filters will be output.
# MIN_SVLEN and MAX_SVLEN can be the string inf or infinity.
#
# The output is a four column tab-delimited table of the form:
# VARIANT_ID SVTYPE SVLEN SAMPLE_ID
# One line will be printed for each sample that has the variant (GT ~ /1/).

BEGIN {
    FS = "\t"
    OFS = "\t"
    CHROM = 1
    POS = 2
    ID = 3
    REF = 4
    ALT = 5
    QUAL = 6
    FILTER = 7
    INFO = 8
    FORMAT = 9
    sfcnt = 0
}

ARGC == 3 && FILENAME == ARGV[1] {
    n = split($0, parts, /;/)
    if (n != 3) {
        print "invalid SV filter: " $0 > "/dev/stderr"
        print "format is SVTYPE;MIN_SVLEN;MAX_SVLEN" > "/dev/stderr"
        exit 84
    }

    if (tolower(parts[2]) ~ /^-?inf(inity)?$/) {
        min = log(0)
    } else {
        min = int(parts[2])
    }

    if (tolower(parts[3]) ~ /^?inf(inity)?$/) {
        max = -log(0)
    } else {
        max = int(parts[3])
    }

    if (min > max) {
        print "invalid SV filter: " $0 > "/dev/stderr"
        print "min SV length cannot be greater than max SV length" > "/dev/stderr"
        exit 85
    }

    sftype[++sfcnt] = parts[1]
    sfmin[sfcnt] = parts[2]
    sfmax[sfcnt] = parts[3]
}

/^#CHROM/ {
    for (i = FORMAT + 1; i < NF; ++i) {
        samples[i] = $i
    }

    next
}

!/^#/ {
    vid = $ID
    q = match($INFO, /;?SVTYPE=[^;]+;?/) 
    svtype = q ? substr($INFO, RSTART, RLENGTH) : ""
    sub(/^;?SVTYPE=/, "", svtype)
    sub(/;$/, "", svtype)
    q = match($INFO, /;?SVLEN=[^;]+;?/)
    svlen = q ? substr($INFO, RSTART, RLENGTH) : ""
    sub(/^;?SVLEN=/, "", svlen)
    sub(/;$/, "", svlen)

    if (svtype == "BND") {
        next
    }

    if (!sfcnt) {
        print_site()
        next
    }

    for (i = 1; i < sfcnt; ++i) {
        if (svtype == sftype[i] && svlen >= sfmin[i] && svlen <= sfmax[i]) {
            print_site()
        }
    }
}

function print_site(    i, gt_fields) {
    for (i = FORMAT + 1; i < NF; ++i) {
        split($i, gt_fields, /:/)
        if (gt_fields[1] ~ /1/) {
            print vid, svtype, svlen, samples[i]
        }
    }
}
