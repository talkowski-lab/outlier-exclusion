# Add OUTLIER filter to FilterGenotypes VCF
# usage: gawk -f flag_outliers.awk <variants> <vcf>
# <variants>
#   List of variant IDs to flag (from JoinRawCalls)
# <vcf>
#   FilterGenotypes VCF (must be streamed if compressed)

BEGIN {
	FS = "\t"
	OFS = "\t"
}

ARGIND == 1 {
	variants[$1]
	next
}

/^##fileformat/ {
	vcf_version = $0
	next
}

/^##contig/ {
	contigs[++ncontigs] = $0
	next
}

/^##/ {
	other_headers[++nheaders] = $0
	next
}

/^##FILTER=<ID=OUTLIER,/ {
	next
}

/^#CHROM/ {
	other_headers[++nheaders] = "##FILTER=<ID=OUTLIER,Description=\"Variant enriched by outlier samples\">"
	nheaders = asort(other_headers)
	print vcf_version
	for (i = 1; i <= nheaders; ++i) {
		print other_headers[i]
	}

	for (i = 1; i <= ncontigs; ++i) {
		print contigs[i]
	}
	# Print the CHROM line
	print $0
	next
}

# The TRUTH_VID value in the FilterGenotypes VCF corresponds to the variant ID
# in the JoinRawCalls VCF.
$8 ~ /TRUTH_VID=/ {
	match($8, /TRUTH_VID=([^;]+);?/, a)
	if (RSTART && a[1] && (a[1] in variants)) {
		if ($7 == "PASS" || $7 == ".") {
			$7 = "OUTLIER"
		} else if ($7 !~ /(^OUTLIER$)|(^OUTLIER;)|(;OUTLIER;)|(;OUTLIER$)/) {
			# Check for the flag first to avoid adding a duplicate.
			# These cases are needed to ensure the flag is not a
			# substring of another flag.
			$7 = $7 ";OUTLIER"
		}
	}
}

1
