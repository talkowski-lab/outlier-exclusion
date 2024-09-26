BEGIN {
    FS = "\t"
}

NR == FNR {
    keys[$1] = $2
    next
}

/\{\{[^{}]+\}\}/ {
    pre = ""
    suf = $0
    match(suf, /\{\{[^{}]+\}\}/)
    while (RSTART) {
        tstart = RSTART
        tend = tstart + RLENGTH - 1
        key = substr(suf, RSTART, RLENGTH)
        gsub(/[{}]/, "", key)

        pre = pre substr(suf, 1, tstart - 1) ((key in keys) ? keys[key] : "")
        suf = substr(suf, tend + 1)
        
        match(suf, /\{\{[^{}]+\}\}/)
    }
    
    print pre suf
    next
}

1