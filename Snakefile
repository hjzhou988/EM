INPUT_DIR = config["input_dir"] # contains {race}.pull and {race}.glid files
OUTPUT_DIR = config["output_dir"]
POPULATIONS = config["populations"]
BLOCKS = ["XR","CB","ACB","XRQ","ACBXRQ"]
Association_file = "config/DRBX_DRB1_assoc.cfg"
JAR = "scripts/em-tools-3.3-2015-08-04.jar"


rule all:
    input:
        # expand(OUTPUT_DIR + "/{block}/{race}.pull", race=POPULATIONS,block=BLOCKS),
        # expand(OUTPUT_DIR + "/{block}/{race}.glid", race=POPULATIONS,block=BLOCKS),
        #expand(OUTPUT_DIR + "/{block}/{race}.nemo", race=POPULATIONS,block=BLOCKS),
        #expand(OUTPUT_DIR + "/{block}/{race}.freqs", race=POPULATIONS,block=BLOCKS ),
        expand(OUTPUT_DIR + "/{block}/{race}.impute",race=POPULATIONS,block=BLOCKS)

rule make_block_2locus:
    input:
        pull_files = INPUT_DIR + "/{race}.pull",
        glid_files = INPUT_DIR + "/{race}.glid"
    output:
        new_pull_files = OUTPUT_DIR + "/{block}/{race}.pull",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
        EM_input_file = OUTPUT_DIR + "/{block}/{race}.nemo"
    wildcard_constraints:
        block="[A-Z]{2}" #Only allow 2-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/make_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/make_block.{block}.{race}.log"
    shell:'''
    perl scripts/make_block-final.pl \
    -i {INPUT_DIR} \
    -o {OUTPUT_DIR} \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -f {Association_file} \
    &> {log}
'''

rule em_block_2locus:
    input:
        EM_input_file = OUTPUT_DIR + "/{block}/{race}.nemo",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
    output:
        freq_file = OUTPUT_DIR + "/{block}/{race}.freqs",
    wildcard_constraints:
        block="[A-Z]{2}" #Only allow 2-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/em_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/em_block.{block}.{race}.log"
    shell:'''
    java -jar -Xmx100g -Dthreads=1 {JAR} w2emo \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -n {input.EM_input_file} \
    -g {input.new_glid_files} \
    -o {output} \
    &> {log}
    '''

rule impute_block_2locus:
    input:
        new_pull_files = OUTPUT_DIR + "/{block}/{race}.pull",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
        freq_file = OUTPUT_DIR + "/{block}/{race}.freqs",
    output:
        OUTPUT_DIR + "/{block}/{race}.impute"
    wildcard_constraints:
        block="[A-Z]{2}" #Only allow 2-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/impute_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/impute_block.{block}.{race}.log"
    shell:'''
    java -jar -Xmx100g -Dthreads=1 {JAR} impute \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -f {input.new_pull_files} \
    -g {input.new_glid_files} \
    -q {input.freq_file} \
    -o {output} \
    &> {log}
    '''

rule make_block_3locus:
    input:
        pull_files = INPUT_DIR + "/{race}.pull",
        glid_files = INPUT_DIR + "/{race}.glid",
        impute_files = expand(OUTPUT_DIR + "/{two_locus_block}/{{race}}.impute",two_locus_block=["CB","XR"])
    output:
        new_pull_files = OUTPUT_DIR + "/{block}/{race}.pull",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
        EM_input_file = OUTPUT_DIR + "/{block}/{race}.nemo"
    wildcard_constraints:
        block="[A-Z]{3}", #Only allow 3-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/make_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/make_block.{block}.{race}.log"
    shell:'''
    perl scripts/make_block-final.pl \
    -i {INPUT_DIR} \
    -o {OUTPUT_DIR} \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -f {Association_file} \
    &> {log}
'''

rule em_block_3locus:
    input:
        EM_input_file = OUTPUT_DIR + "/{block}/{race}.nemo",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
    output:
        freq_file = OUTPUT_DIR + "/{block}/{race}.freqs",
    wildcard_constraints:
        block="[A-Z]{3}" #Only allow 3-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/em_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/em_block.{block}.{race}.log"
    shell:'''
    java -jar -Xmx100g -Dthreads=1 {JAR} w2emo \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -n {input.EM_input_file} \
    -g {input.new_glid_files} \
    -o {output} \
    &> {log}
    '''

rule impute_block_3locus:
    input:
        new_pull_files = OUTPUT_DIR + "/{block}/{race}.pull",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
        freq_file = OUTPUT_DIR + "/{block}/{race}.freqs",
    output:
        OUTPUT_DIR + "/{block}/{race}.impute"
    wildcard_constraints:
        block="[A-Z]{3}" #Only allow 3-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/impute_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/impute_block.{block}.{race}.log"
    shell:'''
    java -jar -Xmx100g -Dthreads=1 {JAR} impute \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -f {input.new_pull_files} \
    -g {input.new_glid_files} \
    -q {input.freq_file} \
    -o {output} \
    &> {log}
    '''
rule make_block_ACBXRQ:
    input:
        pull_files = expand(OUTPUT_DIR + "/{three_locus_block}/{{race}}.pull",three_locus_block=["XRQ"]),
        impute_files = expand(OUTPUT_DIR + "/{three_locus_block}/{{race}}.impute",three_locus_block=["ACB","XRQ"])
    output:
        new_pull_files = OUTPUT_DIR + "/{block}/{race}.pull",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
        EM_input_file = OUTPUT_DIR + "/{block}/{race}.nemo"
    wildcard_constraints:
        block="[A-Z]{6}", #Only allow 6-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/make_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/make_block.{block}.{race}.log"
    shell:'''
    perl scripts/make_block-final.pl \
    -i {INPUT_DIR} \
    -o {OUTPUT_DIR} \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -f {Association_file} \
    &> {log}
'''

rule em_block_ACBXRQ:
    input:
        EM_input_file = OUTPUT_DIR + "/{block}/{race}.nemo",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
    output:
        freq_file = OUTPUT_DIR + "/{block}/{race}.freqs",
    wildcard_constraints:
        block="[A-Z]{6}" #Only allow 6-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/em_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/em_block.{block}.{race}.log"
    shell:'''
    java -jar -Xmx100g -Dthreads=1 {JAR} w2emo \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -n {input.EM_input_file} \
    -g {input.new_glid_files} \
    -o {output} \
    &> {log}
    '''
rule impute_block_ACBXRQ:
    input:
        new_pull_files = OUTPUT_DIR + "/{block}/{race}.pull",
        new_glid_files = OUTPUT_DIR + "/{block}/{race}.glid",
        freq_file = OUTPUT_DIR + "/{block}/{race}.freqs",
    output:
        OUTPUT_DIR + "/{block}/{race}.impute"
    wildcard_constraints:
        block="[A-Z]{6}" #Only allow 6-locus blocks
    benchmark:
        OUTPUT_DIR + "/benchmarks/impute_block.{block}.{race}.benchmark.txt"
    log:
        OUTPUT_DIR + "/logs/impute_block.{block}.{race}.log"
    shell:'''
    java -jar -Xmx100g -Dthreads=1 {JAR} impute \
    -c {wildcards.block} \
    -r {wildcards.race} \
    -f {input.new_pull_files} \
    -g {input.new_glid_files} \
    -q {input.freq_file} \
    -o {output} \
    &> {log}
    '''
