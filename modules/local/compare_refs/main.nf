process COMPARE_REFS {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(ref1), path(ref2)

    output:
    tuple val(meta), path("sha256sum.txt"), emit: txt
    path "versions.yml",                    emit: versions


    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    # Handle gzipped or unzipped file as input.
    # Resolve input to a single uncompressed file in the work dir, using
    # a label prefix to avoid filename collisions (e.g. ref.fa.gz vs ref.fa).
    resolve_gz() {
        local in="\$1"
        local label="\$2"
        local base="\$(basename "\$in")"

        case "\$base" in
            *.gz)
                local out="\${label}_\${base%.gz}"
                gunzip -c "\$in" > "\$out"
                ;;
            *)
                local out="\${label}_\$base"
                cp -L "\$in" "\$out"
                ;;
        esac

        echo "\$out"
    }

    f1="\$(resolve_gz "${ref1}" ref1)"
    f2="\$(resolve_gz "${ref2}" ref2)"

    # Remove ref1_ and ref2_ prefixes from output
    sha256sum "\$f1" "\$f2" | sed 's/\\bref[12]_//' > sha256sum.txt

    # version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sha256sum: "\$(sha256sum --version | sed -n '1p' | cut -f 4 -d ' ')"
        gzip: "\$(gzip --version | sed -n '1p' | cut -f 2 -d ' ')"
    END_VERSIONS
    """
}
