from datetime import date
import pandas as pd
from treetime.utils import numeric_date

segments = ['ha']
lineages = ['h3n2']
path_to_fauna = '../fauna'
resolutions = ['6y', '12y']
passages = ['cell', 'egg']
assays = ['hi', 'fra']
centers = ['who']

wildcard_constraints:
    passage= "(cell|egg)",
    lineage="[a-z0-9]+",
    segment="[a-z]+[0-9]?",
    center="(who)",
    assay="(hi|fra)"


min_length = 900
frequency_regions = ['north_america', 'south_america', 'europe', 'china', 'oceania',
                     'southeast_asia', 'japan_korea', 'south_asia', 'africa']


def vpm(v):
    vpm = {'6m':360, '2y':90, '3y':60, '6y':30, '12y':15}
    return vpm[v.resolution] if v.resolution in vpm else 5

def reference_strain(v):
    references = {'h3n2':"A/Beijing/32/1992",
                  'h1n1pdm':"A/California/07/2009",
                  'vic':"B/HongKong/02/1993",
                  'yam':"B/Singapore/11/1994"
                  }
    return references[v.lineage]

genes_to_translate = {'ha':['SigPep', 'HA1', 'HA2'], 'na':['NA']}
def gene_names(w):
    return genes_to_translate[w.segment]

def translations(w):
    genes = gene_names(w)
    return ["results/aa-seq_%s_%s_%s_%s_%s_%s_%s.fasta"%(w.center, w.lineage, w.segment, w.resolution, w.passage, w.assay, g)
            for g in genes]

def pivot_interval(w):
    """Returns the number of months between pivots by build resolution.
    """
    pivot_intervals_by_resolution = {'6m': 1, '2y': 1, '3y': 2, '6y': 3, '12y': 6}
    return pivot_intervals_by_resolution[w.resolution]

def min_date(wildcards):
    now = numeric_date(date.today())
    if wildcards.resolution[-1] == "y":
        years_back = int(wildcards.resolution[:-1])
    elif wildcards.resolution[-1] == "m":
        years_back = int(wildcards.resolution[:-1]) / 12.
    else:
        years_back = 3
    return now - years_back

def max_date(w):
    return numeric_date(date.today())

def clock_rate(w):
    rate = {
        ('h3n2', 'ha'): 0.0043, ('h3n2', 'na'):0.0029,
        ('h1n1pdm', 'ha'): 0.0040, ('h1n1pdm', 'na'):0.0032,
        ('vic', 'ha'): 0.0024, ('vic', 'na'):0.0015,
        ('yam', 'ha'): 0.0019, ('yam', 'na'):0.0013
    }
    return rate[(w.lineage, w.segment)]

#
# Define clades functions
#
def _get_clades_file_for_wildcards(wildcards):
    if wildcards.segment == "ha":
        return "config/clades_%s_ha.tsv"%(wildcards.lineage)
    else:
        return "results/clades_%s_%s_ha_%s_%s_%s.json"%(wildcards.center, wildcards.lineage,
                                                        wildcards.resolution, wildcards.passage, wildcards.assay)


#
# Define titer data sets to be used.
#
def _get_tdb_databases(wildcards):
    if wildcards.center == 'who':
        return "cdc_tdb crick_tdb niid_tdb vidrl_tdb tdb"
    else:
        return wildcards.center + "_tdb tdb"

def exclude_where(wildcards):
    return "country=? region=?"


#
# Define LBI parameters and functions.
#
LBI_params = {
    '6m': {"tau": 0.3, "time_window": 0.5},
    '2y': {"tau": 0.3, "time_window": 0.5},
    '3y': {"tau": 0.4, "time_window": 0.6},
    '6y': {"tau": 0.25, "time_window": 0.75},
    '12y': {"tau": 0.25, "time_window": 0.75}
}

def _get_lbi_tau_for_wildcards(wildcards):
    return LBI_params[wildcards.resolution]["tau"]

def _get_lbi_window_for_wildcards(wildcards):
    return LBI_params[wildcards.resolution]["time_window"]

#
# Configure amino acid distance masks.
#

# Load mask configuration including which masks map to which attributes per
# lineage and segment.
masks_config = pd.read_table("config/mask_config.tsv")

def _get_build_mask_config(wildcards):
    config = masks_config[(masks_config["lineage"] == wildcards.lineage) &
                          (masks_config["segment"] == wildcards.segment)]
    if config.shape[0] > 0:
        return config
    else:
        return None

def _get_mask_attribute_names_by_wildcards(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "attribute"].values)

def _get_mask_names_by_wildcards(wildcards):
    config = _get_build_mask_config(wildcards)
    return " ".join(config.loc[:, "mask"].values)


rule all_egg_analyses:
    input:
        dataframe = expand("dataframes/{lineage}_{segment}_{resolution}_{assay}.csv", lineage=lineages, segment=segments, resolution=resolutions, assay=assays),
        aa_mut_plot = expand("plots/{lineage}_{segment}_{resolution}_{assay}/egg_mutation_aa_prevalence_{lineage}_{segment}_{resolution}_{assay}.pdf", lineage=lineages, segment=segments, resolution=resolutions, assay=assays),
        epistasis_plot = expand("plots/{lineage}_{segment}_{resolution}_{assay}/epistasis_heatmap_{lineage}_{segment}_{resolution}_{assay}.pdf", lineage=lineages, segment=segments, resolution=resolutions, assay=assays),
        chord_plot = expand("plots/{lineage}_{segment}_{resolution}_{assay}/epistasis_chord_diagram_{lineage}_{segment}_{resolution}_{assay}.pdf", lineage=lineages, segment=segments, resolution=resolutions, assay=assays),
        pairs_json = expand("egg_results/egg_mutation_accuracy_{lineage}_{segment}_{resolution}_{assay}.json", lineage=lineages, segment=segments, resolution=resolutions, assay=assays)

rule all:
    input:
        auspice_tree = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tree.json",
                              lineage=lineages, segment=segments, resolution=resolutions,
                              assay=assays),
        auspice_meta = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_meta.json",
                              lineage=lineages, segment=segments, resolution=resolutions,
                              assay=assays),
        auspice_tip_frequencies = expand("auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tip-frequencies.json",
                              lineage=lineages, segment=segments, resolution=resolutions,
                              assay=assays)

rule download_all:
    input:
        titers = expand("data/{center}_{lineage}_{passage}_{assay}_titers.tsv",
                         center=centers, lineage=lineages, passage=passages, assay=assays),
        sequences = expand("data/{lineage}_{segment}.fasta", lineage=lineages, segment=segments),
        concatenated_titers = expand("data/{center}_{lineage}_concat_{assay}_titers.tsv",
                                    center=centers, lineage=lineages, assay=assays),
        metadata = expand("results/metadata_{lineage}_{segment}.tsv", lineage=lineages,
                            segment=segments)

def _get_node_data_for_export(wildcards):
    """Return a list of node data files to include for a given build's wildcards.
    """
    # Define inputs shared by all builds.
    inputs = [
        rules.refine.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.titers_tree.output.titers_model,
        rules.titers_sub.output.titers_model,
        rules.clades.output.clades,
        rules.traits.output.node_data,
        rules.lbi.output.lbi
    ]

    # Only request a distance file for builds that have mask configurations
    # defined.
    if _get_build_mask_config(wildcards) is not None:
        inputs.append(rules.distances.output.distances)

    # Convert input files from wildcard strings to real file names.
    inputs = [input_file.format(**wildcards) for input_file in inputs]
    return inputs

rule files:
    params:
        outliers = "config/outliers_{lineage}.txt",
        exclude_sites = "config/exclude-sites_{lineage}.txt",
        references = "config/references_{lineage}.txt",
        force_include = "config/forceinclude_{lineage}.txt",
        reference = "config/reference_{lineage}_{segment}.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config_{lineage}.json",

files = rules.files.params


rule download_sequences:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{lineage}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus flu \
            --fasta_fields {params.fasta_fields} \
            --resolve_method split_passage \
            --select locus:{wildcards.segment} lineage:seasonal_{wildcards.lineage} \
            --path data \
            --fstem {wildcards.lineage}_{wildcards.segment} \

        mv data/{wildcards.lineage}_{wildcards.segment}.fasta  \
            data/{wildcards.lineage}_{wildcards.segment}_temp.fasta \

        python3 scripts/tag_cell_seqs.py \
            --in_file data/{wildcards.lineage}_{wildcards.segment}_temp.fasta \
            --out_file {output.sequences} \
        """

rule download_titers:
    message: "Downloading titers from fauna: {wildcards.lineage}, {wildcards.assay}, {wildcards.center}"
    output:
        titers = "data/{center}_{lineage}_{passage}_{assay}_titers.tsv"
    params:
        dbs = _get_tdb_databases,
    shell:
        """
        python3 {path_to_fauna}/tdb/download.py \
            --database {params.dbs} \
            --virus flu \
            --subtype {wildcards.lineage} \
            --select assay_type:{wildcards.assay} serum_passage_category:{wildcards.passage} \
            --path data \
            --fstem {wildcards.center}_{wildcards.lineage}_{wildcards.passage}_{wildcards.assay}
        """

rule concatenate_titers:
    message: "Concatenating egg and cell titers"
    input:
        titer_tsvs = expand("data/{{center}}_{{lineage}}_{passage}_{{assay}}_titers.tsv", passage=passages)
    output:
        concatenated_titers = "data/{center}_{lineage}_concat_{assay}_titers.tsv",
    shell:
        """
        python3 scripts/concatenate_titers.py \
            --titer_tsvs {input.titer_tsvs} \
            --output {output.concatenated_titers}
        """

passages = ['concat']
wildcard_constraints:
    passage= "(concat)",
    lineage="[a-z0-9]+",
    segment="[a-z]+[0-9]?",
    center="(who)",
    assay="(hi|fra)"

rule find_egg_seqs:
    message: "Finding all egg-passaged viruses"
    input:
        sequences = expand("data/{lineage}_{segment}.fasta", lineage=lineages, segment=segments),
        references = files.references
    output:
        force_include = "config/forceinclude_{lineage}.txt"
    shell:
        """
        python3 scripts/find_egg_seqs.py \
            --fasta_file {input.sequences} \
            --references_file {input.references} \
            --force_include_file {output.force_include} \
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download_sequences.output.sequences
    output:
        sequences = "results/sequences_{lineage}_{segment}.fasta",
        metadata = "results/metadata_{lineage}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus accession date region country division location passage authors age gender"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message:
        """
        Filtering {wildcards.lineage} {wildcards.segment} {wildcards.passage} sequences:
          - less than {params.min_length} bases
          - outliers
          - samples with missing region and country metadata
          - samples that are egg-passaged if cell build
        """
    input:
        metadata = rules.parse.output.metadata,
        sequences = rules.parse.output.sequences,
        exclude = files.outliers
    output:
        sequences = 'results/filtered_{lineage}_{segment}_{passage}.fasta'
    params:
        min_length = min_length,
        exclude_where = exclude_where
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --min-length {params.min_length} \
            --non-nucleotide \
            --exclude {input.exclude} \
            --exclude-where {params.exclude_where} \
            --output {output}
        """

rule select_strains:
    input:
        sequences = expand("results/filtered_{{lineage}}_{segment}_{{passage}}.fasta", segment=segments),
        metadata = expand("results/metadata_{{lineage}}_{segment}.tsv", segment=segments),
        titers = rules.concatenate_titers.output.concatenated_titers,
        include = files.force_include
    output:
        strains = "results/strains_{center}_{lineage}_{resolution}_{passage}_{assay}.txt",
    params:
        viruses_per_month = vpm
    shell:
        """
        python3 scripts/select_strains.py \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --segments {segments} \
            --include {input.include} \
            --lineage {wildcards.lineage} \
            --resolution {wildcards.resolution} \
            --viruses_per_month {params.viruses_per_month} \
            --titers {input.titers} \
            --output {output.strains}
        """

rule extract:
    input:
        sequences = rules.filter.output.sequences,
        strains = rules.select_strains.output.strains
    output:
        sequences = 'results/extracted_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.fasta'
    shell:
        """
        python3 scripts/extract_sequences.py \
            --sequences {input.sequences} \
            --samples {input.strains} \
            --output {output}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
        """
    input:
        sequences = rules.extract.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --nthreads 1
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment,
        exclude_sites = files.exclude_sites
    output:
        tree = "results/tree-raw_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads 1 \
            --exclude-sites {input.exclude_sites}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.nwk",
        node_data = "results/branch-lengths_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock_rate = clock_rate
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule reconstruct_translations:
    message: "Reconstructing translations required for titer models and frequencies"
    input:
        tree = rules.refine.output.tree,
        node_data = "results/aa-muts_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    output:
        aa_alignment = "results/aa-seq_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_{gene}.fasta"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene {wildcards.gene} \
            --output {output.aa_alignment} \
            --internal-nodes
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    params:
        columns = "region"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule titers_sub:
    input:
        titers = rules.concatenate_titers.output.concatenated_titers,
        aa_muts = rules.translate.output,
        alignments = translations,
        tree = rules.refine.output.tree
    params:
        genes = gene_names
    output:
        titers_model = "results/titers-sub-model_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    shell:
        """
        augur titers sub \
            --titers {input.titers} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule titers_tree:
    input:
        titers = rules.concatenate_titers.output.concatenated_titers,
        tree = rules.refine.output.tree
    output:
        titers_model = "results/titers-tree-model_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    shell:
        """
        augur titers tree \
            --titers {input.titers} \
            --tree {input.tree} \
            --output {output.titers_model}
        """

rule mutation_frequencies:
    input:
        metadata = rules.parse.output.metadata,
        alignment = translations
    params:
        genes = gene_names,
        min_date = min_date,
        max_date = max_date,
        pivot_interval = pivot_interval
    output:
        mut_freq = "results/mutation-frequencies_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    shell:
        """
        augur frequencies \
            --alignments {input.alignment} \
            --metadata {input.metadata} \
            --gene-names {params.genes} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --pivot-interval {params.pivot_interval} \
            --output {output.mut_freq}
        """

rule tip_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        weights = "config/frequency_weights_by_region.json"
    params:
        narrow_bandwidth = 1 / 12.0,
        wide_bandwidth = 3 / 12.0,
        proportion_wide = 0.0,
        weight_attribute = "region",
        min_date = min_date,
        max_date = max_date,
        pivot_interval = pivot_interval
    output:
        tip_freq = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_tip-frequencies.json"
    shell:
        """
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --proportion-wide {params.proportion_wide} \
            --weights {input.weights} \
            --weights-attribute {params.weight_attribute} \
            --pivot-interval {params.pivot_interval} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output}
        """

rule tree_frequencies:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
    params:
        min_date = min_date,
        max_date = max_date,
        pivot_interval = pivot_interval,
        regions = ['global'] + frequency_regions,
        min_clade = 20
    output:
        "results/tree-frequencies_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    shell:
        """
        augur frequencies \
            --method diffusion \
            --include-internal-nodes \
            --tree {input.tree} \
            --regions {params.regions} \
            --metadata {input.metadata} \
            --pivot-interval {params.pivot_interval} \
            --minimal-clade-size {params.min_clade} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --output {output}
        """


rule clades:
    message: "Annotating clades"
    input:
        tree = "results/tree_{center}_{lineage}_ha_{resolution}_{passage}_{assay}.nwk",
        nt_muts = rules.ancestral.output,
        aa_muts = rules.translate.output,
        clades = _get_clades_file_for_wildcards
    output:
        clades = "results/clades_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    run:
        if wildcards.segment == 'ha':
            shell("""
                augur clades \
                    --tree {input.tree} \
                    --mutations {input.nt_muts} {input.aa_muts} \
                    --clades {input.clades} \
                    --output {output.clades}
            """)
        else:
            shell("""
                python3 scripts/import_tip_clades.py \
                    --tree {input.tree} \
                    --clades {input.clades} \
                    --output {output.clades}
            """)

rule distances:
    input:
        tree = rules.refine.output.tree,
        alignments = translations,
        masks = "config/{segment}_masks.tsv"
    params:
        genes = gene_names,
        attribute_names = _get_mask_attribute_names_by_wildcards,
        mask_names = _get_mask_names_by_wildcards
    output:
        distances = "results/distances_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json",
    shell:
        """
        augur distance \
            --tree {input.tree} \
            --alignment {input.alignments} \
            --gene-names {params.genes} \
            --masks {input.masks} \
            --output {output} \
            --attribute-names {params.attribute_names} \
            --mask-names {params.mask_names}
        """

rule lbi:
    message: "Calculating LBI"
    input:
        tree = rules.refine.output.tree,
        branch_lengths = rules.refine.output.node_data
    params:
        tau = _get_lbi_tau_for_wildcards,
        window = _get_lbi_window_for_wildcards,
        names = "lbi"
    output:
        lbi = "results/lbi_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}.json"
    shell:
        """
        augur lbi \
            --tree {input.tree} \
            --branch-lengths {input.branch_lengths} \
            --output {output} \
            --attribute-names {params.names} \
            --tau {params.tau} \
            --window {params.window}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        auspice_config = files.auspice_config,
        node_data = _get_node_data_for_export
    output:
        auspice_tree = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_tree.json",
        auspice_meta = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_meta.json",
        auspice_seq = "auspice/flu_{center}_{lineage}_{segment}_{resolution}_{passage}_{assay}_root-sequence.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta} \
            --output-sequence {output.auspice_seq}
        """

rule simplify_auspice_names:
    input:
        tree = "auspice/flu_who_{lineage}_{segment}_{resolution}_concat_{assay}_tree.json",
        meta = "auspice/flu_who_{lineage}_{segment}_{resolution}_concat_{assay}_meta.json",
        seq = "auspice/flu_who_{lineage}_{segment}_{resolution}_concat_{assay}_root-sequence.json",
        frequencies = "auspice/flu_who_{lineage}_{segment}_{resolution}_concat_{assay}_tip-frequencies.json"
    output:
        tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tree.json",
        meta = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_meta.json",
        seq = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_root-sequence.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tip-frequencies.json"
    shell:
        '''
        mv {input.tree} {output.tree} &
        mv {input.meta} {output.meta} &
        mv {input.seq} {output.seq} &
        mv {input.frequencies} {output.frequencies} &
        '''

rule organize_output:
    input:
        tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tree.json",
        seqs = "results/aa-seq_who_{lineage}_{segment}_{resolution}_concat_{assay}_HA1.fasta",
        root_seq = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_root-sequence.json",
    output:
        dataframe = "dataframes/{lineage}_{segment}_{resolution}_{assay}.csv",
        egg_dataframe = "dataframes/{lineage}_{segment}_{resolution}_{assay}_egg.csv",
        tidy_dataframe = "dataframes/{lineage}_{segment}_{resolution}_{assay}_tidy.csv",

    shell:
        """
        python3 scripts/organize_output.py \
            --tree {input.tree} \
            --seqs {input.seqs} \
            --root_seq {input.root_seq} \
        """

rule find_egg_mutations:
    input:
        dataframe = "dataframes/{lineage}_{segment}_{resolution}_{assay}.csv",
    output:
        aa_mut_plot = "plots/{lineage}_{segment}_{resolution}_{assay}/egg_mutation_aa_prevalence_{lineage}_{segment}_{resolution}_{assay}.pdf",
        site_mut_plot = "plots/{lineage}_{segment}_{resolution}_{assay}/egg_mutation_site_prevalence_{lineage}_{segment}_{resolution}_{assay}.pdf",
        before_after_plot = "plots/{lineage}_{segment}_{resolution}_{assay}/before_after_eggpassaging_{lineage}_{segment}_{resolution}_{assay}.pdf",
    shell:
        """
        python3 scripts/find_egg_mutations.py \
            --in_file {input.dataframe} \
        """

rule egg_epistasis:
    input:
        dataframe = "dataframes/{lineage}_{segment}_{resolution}_{assay}.csv",
    output:
        epistasis_plot = "plots/{lineage}_{segment}_{resolution}_{assay}/epistasis_heatmap_{lineage}_{segment}_{resolution}_{assay}.pdf",
        chord_plot = "plots/{lineage}_{segment}_{resolution}_{assay}/epistasis_chord_diagram_{lineage}_{segment}_{resolution}_{assay}.pdf",
    shell:
        """
        python3 scripts/plot_egg_epistasis.py \
            --in_file {input.dataframe} \
        """

rule compair_pairs:
    input:
        dataframe = "dataframes/{lineage}_{segment}_{resolution}_{assay}.csv",
    output:
        pairs_json = "egg_results/egg_mutation_accuracy_{lineage}_{segment}_{resolution}_{assay}.json",
    shell:
        """
        python3 scripts/compare_pairs.py \
            --in_file {input.dataframe} \
        """

rule targets:
    input:
        tree = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tree.json",
        meta = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_meta.json",
        seq = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_root-sequence.json",
        frequencies = "auspice/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}_tip-frequencies.json"
    output:
        target = "targets/flu_seasonal_{lineage}_{segment}_{resolution}_{assay}"
    shell:
        '''
        touch {output.target}
        '''

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "targets ",
        "auspice ",
        "auspice-who"
    shell:
        "rm -rfv {params}"
