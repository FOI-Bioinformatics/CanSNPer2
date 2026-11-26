#!/usr/bin/env python3
"""
cansnper2 - one command with a verb for every job (2.6.0).

Verbs: cansnper2 type, selftest, db create|import|export|info|convert|modify|
remove, scheme list|get, ref get. They call the run functions of the commands
before 2.6.0 (CanSNPer2, CanSNPer2-database, CanSNPer2-download, CanSNPer2-test),
which are no longer installed. The names after every dest= are the internal
names those run functions read.
"""

import argparse
import logging

from CanSNPer2 import __version__
from CanSNPer2.modules.TopLevel import run_main
from CanSNPer2.modules.Aligners import ALIGNERS, DEFAULT_ALIGNER
from CanSNPer2.modules.ParsePileup import PLATFORMS, DEFAULT_PLATFORM

import CanSNPer2.CanSNPerTree as typing_command
import CanSNPer2.SNPDatabase as database_command
import CanSNPer2.DownloadGenomes as download_command
import CanSNPer2.selftest as selftest_command


class HelpFormatter(argparse.HelpFormatter):
    """help that appends "(default: X)" only to options that take a value and have a
    default worth showing. Flags, options without a default and help texts that
    already state their default are left as written; ArgumentDefaultsHelpFormatter
    printed "(default: None)", "(default: False)" and the internal value of --quiet.
    A description or epilog with line breaks (the examples) is printed as written,
    one without them is wrapped."""

    def _fill_text(self, text, width, indent):
        if "\n" in text:
            return "".join(indent + line for line in text.splitlines(keepends=True))
        return super()._fill_text(text, width, indent)

    def _get_help_string(self, action):
        text = action.help or ""
        if "default" in text or not action.option_strings or action.nargs == 0:
            return text
        ## a 0 default is shown; 0 == False, so False is tested by identity
        if action.default is None or action.default is False or action.default in ("", [], argparse.SUPPRESS):
            return text
        return text + " (default: %(default)s)"


def examples(*lines):
    """the epilog of a verb: a few command lines"""
    return "examples:\n" + "\n".join("  " + line for line in lines)


def add_logging_options(parser, logs_dest, logs_default, logs_help="folder for the log files"):
    """-v/--debug/--quiet and the log folder, the same three levels the old commands
    have (--supress becomes --quiet); logs_dest and logs_default are the internal
    name and default of the log folder of this command"""
    group = parser.add_argument_group("Logging options")
    group.add_argument("--log-dir", dest=logs_dest, metavar="DIR", default=logs_default, help=logs_help)
    group.add_argument(
        "-v", "--verbose", action="store_const", const=logging.INFO, help="print progress information"
    )
    group.add_argument(
        "--debug",
        action="store_const",
        const=logging.DEBUG,
        help="print debug information, and the traceback of an unexpected error",
    )
    group.add_argument(
        "--quiet",
        dest="supress",
        action="store_const",
        const=logging.ERROR,
        default=logging.WARNING,
        help="print errors only, not warnings",
    )


def add_scheme_source_options(parser):
    """where schemes are read from: CanSNPer2-data on GitHub, another repository,
    or a local folder laid out the same way"""
    group = parser.add_argument_group("Scheme source options")
    group.add_argument(
        "--repo",
        metavar="OWNER/NAME",
        help="GitHub repository laid out as CanSNPer2-data (default: FOI-Bioinformatics/CanSNPer2-data)",
    )
    group.add_argument(
        "--branch",
        metavar="REF",
        help="branch, tag or commit of the repository; in CanSNPer2-data, master "
        "holds the verified schemes and dev also the draft schemes (default: "
        "master, or the default branch of the repository given with --repo)",
    )
    group.add_argument(
        "--data-dir",
        dest="data_dir",
        metavar="DIR",
        help="local directory laid out as CanSNPer2-data to read the schemes from "
        "instead of GitHub; cannot be combined with --repo or --branch",
    )


def add_type_subparser(subparsers):
    parser = subparsers.add_parser(
        "type",
        help="type assemblies or sequencing reads against a scheme database",
        description="Type assemblies or sequencing reads against a scheme database: "
        "each query is compared with the reference genomes of the database at the "
        "canSNP positions, and the deepest node of the scheme tree that the derived "
        "canSNPs support is reported.",
        epilog=examples(
            "cansnper2 type sample1.fna sample2.fna -d francisella_tularensis.db",
            "cansnper2 type queries.txt -d scheme.db --ref-dir references/ -o results/",
            "cansnper2 type s1_R1.fastq.gz s1_R2.fastq.gz -d scheme.db \\",
            "    --platform illumina",
            "cansnper2 type --samplesheet samples.tsv -d scheme.db",
        ),
        formatter_class=HelpFormatter,
    )
    required = parser.add_argument_group("Input")
    required.add_argument(
        "query",
        nargs=argparse.ZERO_OR_MORE,
        metavar="query",
        help="assemblies (FASTA) or sequencing reads (FASTQ, may be gzip-compressed) "
        "to type, or one list file with one file path per line; paired read "
        "files are recognized by their names (_R1/_R2, _1/_2)",
    )
    required.add_argument("-d", "--database", metavar="DB", help="the scheme database (.db); required")

    output = parser.add_argument_group("Output options")
    output.add_argument(
        "-o",
        "--outdir",
        metavar="DIR",
        help="output directory; the overview summary.html is written there and "
        "the per-sample files in its samples folder",
    )
    output.add_argument(
        "--save-tree",
        dest="save_tree",
        action="store_true",
        help="also write a standalone tree page per sample "
        "(samples/<sample>_tree.html, the overview already shows every tree)",
    )
    output.add_argument(
        "--no-snp-files",
        dest="no_snpfiles",
        action="store_false",
        help="do not write the per-sample canSNP files; the call is still printed "
        "and written to the run summary",
    )
    output.add_argument(
        "--summary",
        action="store_true",
        help="also write snps_summary.txt and a tree with the final SNPs of all samples",
    )
    output.add_argument(
        "--node-labels",
        dest="annotations",
        metavar="FILE",
        help="tab-separated file of labels for tree nodes in the HTML pages: "
        "node name, kind (node or clade), rank, label[, color]",
    )

    run = parser.add_argument_group("Run options")
    run.add_argument(
        "--ref-dir",
        dest="refdir",
        metavar="DIR",
        help="folder with the reference genomes of the database (<genome>.fna)",
    )
    run.add_argument(
        "--allow-missing-references",
        dest="allow_missing_references",
        action="store_true",
        help="type even when reference genomes have no file in --ref-dir",
    )
    run.add_argument(
        "--workdir",
        metavar="DIR",
        help="directory that relative paths of the output, log and temporary folders are resolved against",
    )
    run.add_argument(
        "--min-required-hits",
        dest="min_required_hits",
        metavar="N",
        type=typing_command.positive_int,
        help="minimum number of derived canSNPs on the path from the root to a "
        "node for that node to be called",
    )
    run.add_argument(
        "--strictness",
        metavar="FRACTION",
        type=typing_command.fraction,
        help="minimum fraction of the canSNPs on that path that must be derived; "
        "the others may lack data (0-1)",
    )
    run.add_argument(
        "--mask",
        metavar="BASES",
        type=typing_command.non_negative_int,
        help="assemblies only: ignore canSNPs closer than this many bases to the edge of an alignment block",
    )
    run.add_argument(
        "--keep-going",
        dest="keep_going",
        action="store_true",
        help="continue with the other queries when one cannot be typed",
    )
    run.add_argument(
        "--rerun", action="store_true", help="type queries again although their result file exists"
    )
    run.add_argument(
        "--aligner", choices=sorted(ALIGNERS), help="program that aligns an assembly to each reference genome"
    )
    run.add_argument(
        "--alignment-timeout",
        dest="alignment_timeout",
        type=typing_command.minutes,
        metavar="MINUTES",
        help="stop an alignment that runs longer than this; 0 for no limit",
    )
    run.add_argument(
        "--aligner-path",
        dest="aligner_path",
        metavar="PATH",
        help="the aligner program, or the folder that holds it (default: searched in PATH)",
    )
    run.add_argument(
        "--skip-alignment",
        dest="skip_alignment",
        action="store_true",
        help="do not align; read the alignments kept in --tmp-dir by an earlier run with --keep-temp",
    )
    run.add_argument(
        "--keep-temp", dest="keep_temp", action="store_true", help="keep the temporary alignment files"
    )
    run.add_argument(
        "--tmp-dir",
        dest="tmpdir",
        metavar="DIR",
        default=None,
        help="folder for temporary alignments (default: CanSNPer2 in the "
        "system temporary folder given by TMPDIR, or <outdir>/.cansnper2_tmp "
        "when that cannot be written)",
    )

    add_read_options(parser)
    add_logging_options(parser, logs_dest="logdir", logs_default="logs/")
    """the internal options the new interface does not offer keep their old
		defaults: --read_input and the hidden --no_save_tree were dropped,
		--version is a global option of cansnper2"""
    parser.set_defaults(
        read_input=False,
        no_save_tree=False,
        version=False,
        query=[],
        database=None,
        outdir="results",
        save_tree=False,
        no_snpfiles=True,
        summary=False,
        annotations=None,
        refdir="references/",
        allow_missing_references=False,
        workdir="./",
        min_required_hits=3,
        strictness=0.7,
        mask=0,
        keep_going=False,
        rerun=False,
        aligner=DEFAULT_ALIGNER,
        alignment_timeout=60,
        aligner_path="",
        skip_alignment=False,
        keep_temp=False,
    )
    return parser


def add_read_options(parser):
    """read typing and local references: the queries are sequencing reads, which
    minimap2 maps to each reference and samtools mpileup counts at the canSNPs"""
    reads = parser.add_argument_group(
        "Read typing options",
        "used when the queries are sequencing reads (minimap2 and samtools are "
        "required); thresholds that are not given take the defaults of --platform",
    )
    reads.add_argument(
        "--platform",
        choices=sorted(PLATFORMS),
        default=DEFAULT_PLATFORM,
        help="sequencing platform: illumina (short reads, paired or single), "
        "ont (Oxford Nanopore) or hifi (PacBio HiFi)",
    )
    reads.add_argument(
        "--samplesheet",
        metavar="FILE",
        default=None,
        help="tab-separated file with the header columns sample, reads1 and "
        "optionally reads2 and platform, one sample per row; replaces the "
        "query arguments",
    )
    reads.add_argument(
        "--min-depth",
        dest="min_depth",
        metavar="N",
        type=typing_command.positive_int,
        default=None,
        help="minimum number of A, C, G or T bases at a canSNP position, after the "
        "mapping and base quality filters, for its state to be read (default: 5 "
        "for all platforms)",
    )
    reads.add_argument(
        "--min-allele-fraction",
        dest="min_allele_fraction",
        metavar="FRACTION",
        type=typing_command.fraction,
        default=None,
        help="minimum fraction of the reads that must carry one base for the "
        "site to take that state; below it the site is reported as mixed "
        "(default: illumina 0.8, ont 0.7, hifi 0.8)",
    )
    reads.add_argument(
        "--min-mapq",
        dest="min_mapq",
        metavar="Q",
        type=typing_command.non_negative_int,
        default=None,
        help="minimum mapping quality of a read (default: 20 for all platforms)",
    )
    reads.add_argument(
        "--min-baseq",
        dest="min_baseq",
        metavar="Q",
        type=typing_command.non_negative_int,
        default=None,
        help="minimum base quality (default: illumina 20, ont 10, hifi 20)",
    )
    reads.add_argument(
        "--mixed-sites",
        dest="mixed_sites",
        metavar="N",
        type=typing_command.positive_int,
        default=3,
        help="number of mixed canSNP sites from which a sample is flagged as possibly mixed",
    )
    reads.add_argument(
        "--threads",
        metavar="N",
        type=typing_command.positive_int,
        default=4,
        help="threads for minimap2 and samtools sort, per sample; they are divided "
        "between the reference genomes the reads are mapped to at the same time, "
        "at least one each",
    )
    reads.add_argument(
        "--max-reads",
        dest="max_reads",
        metavar="N",
        type=typing_command.non_negative_int,
        default=None,
        help="number of reads mapped per sample, the first ones in each read file; "
        "read pairs for paired-end reads; 0 maps all reads (default: 300000 for "
        "all platforms)",
    )
    reads.add_argument(
        "--local-references",
        dest="local_references",
        action="store_true",
        help="map the reads again to the representative genome of the scheme "
        "closest to the sample and read the canSNPs at their positions in that "
        "genome (genomes in <ref-dir>/local; see ref get --local-references)",
    )
    reads.add_argument(
        "--local-rounds",
        dest="local_rounds",
        metavar="N",
        type=typing_command.positive_int,
        default=2,
        help="with --local-references: maximum number of rounds of choosing a closer representative genome",
    )
    reads.add_argument(
        "--screen-reads",
        dest="screen_reads",
        metavar="N",
        type=typing_command.positive_int,
        default=100000,
        help="with --local-references: number of reads per file that are mapped "
        "to the candidate genomes to choose the local reference",
    )
    reads.add_argument(
        "--fetch-local",
        dest="fetch_local",
        action="store_true",
        help="with --local-references: download from NCBI during typing the "
        "representative genomes that are needed and not yet in <ref-dir>/local",
    )


def database_defaults(parser, **mode):
    """every internal option of the database command with its default, and the
    one mode of this verb set to True"""
    settings = dict(
        database=None,
        tree=False,
        annotation=False,
        references=False,
        source_type="CanSNPer",
        create=False,
        do_import=False,
        scheme_name=None,
        scheme_status=None,
        scheme_source=None,
        info=False,
        convert=False,
        refdir=None,
        mod_file=None,
        parent=None,
        remove=False,
        replace=False,
        export=False,
        export_format=None,
        outdir=None,
        tmpdir="/tmp",
        version=False,
        genomes=None,
    )
    settings.update(mode)
    parser.set_defaults(**settings)


def add_scheme_info_options(parser):
    """the scheme identity recorded in a new database"""
    group = parser.add_argument_group("Scheme information (recorded in the database)")
    group.add_argument(
        "--scheme-name",
        dest="scheme_name",
        metavar="NAME",
        help="name of the scheme, printed by db info and in the typing results",
    )
    group.add_argument(
        "--scheme-status", dest="scheme_status", metavar="STATUS", help="verified, draft or unknown"
    )
    group.add_argument(
        "--scheme-source",
        dest="scheme_source",
        metavar="SRC",
        help="where the scheme files came from, e.g. a repository and commit",
    )


def add_source_file_options(parser, required):
    """the three source files of db create and db import; --annotation becomes
    --cansnps because the 2.4.0 file it reads lists canSNPs. db import needs all
    three, db create builds a database from any of them"""
    needed = "; required" if required else ""
    parser.add_argument("--tree", metavar="FILE", help="tree file of the scheme (<name>_tree.txt)" + needed)
    parser.add_argument(
        "--cansnps",
        dest="annotation",
        metavar="FILE",
        help="canSNP file of the scheme (<name>_canSNPs.txt)" + needed,
    )
    parser.add_argument(
        "--references",
        metavar="FILE",
        help="reference genomes file of the scheme (<name>_references.txt)" + needed,
    )
    parser.add_argument(
        "--genomes",
        metavar="FILE",
        help="genomes file of the scheme (<name>_genomes.txt: genome, node, role); "
        "optional, used by read typing with --local-references",
    )


def add_db_subparsers(subparsers):
    db = subparsers.add_parser(
        "db",
        description="create, import, export, inspect and modify CanSNP databases",
        help="create, import, export, info, convert, modify, remove",
    )
    db_sub = db.add_subparsers(dest="db_action", metavar="<action>", required=True)
    parsers = {}

    create = db_sub.add_parser(
        "create",
        help="build a database from the files of a scheme",
        description="Build a database from the tree, canSNP and reference genome "
        "files of a scheme; a database can also be started from only some of "
        "them and completed with db modify. The files are checked first; the "
        "database is written only when they pass.",
        epilog=examples(
            "cansnper2 db create -d my_scheme.db --tree my_tree.txt \\",
            "    --cansnps my_canSNPs.txt --references my_references.txt \\",
            "    --scheme-name my_scheme --scheme-status draft",
        ),
        formatter_class=HelpFormatter,
    )
    create.add_argument("-d", "--database", metavar="DB", help="the new database to write; required")
    add_source_file_options(create, required=False)
    add_scheme_info_options(create)
    add_logging_options(create, logs_dest="logs", logs_default="logs")
    parsers["db create"] = create

    imp = db_sub.add_parser(
        "import",
        help="build a database from the files written by db export",
        description="Build a database from the files written by db export; the "
        "same checks and options as db create.",
        epilog=examples(
            "cansnper2 db import -d new.db --tree export/scheme_tree.txt \\",
            "    --cansnps export/scheme_canSNPs.txt \\",
            "    --references export/scheme_references.txt",
        ),
        formatter_class=HelpFormatter,
    )
    imp.add_argument("-d", "--database", metavar="DB", help="the new database to write; required")
    add_source_file_options(imp, required=True)
    add_scheme_info_options(imp)
    add_logging_options(imp, logs_dest="logs", logs_default="logs")
    parsers["db import"] = imp

    info = db_sub.add_parser(
        "info",
        help="print the scheme information and size of a database",
        description="Print the scheme information a database records (name, status, "
        "source, the CanSNPer2 version that wrote it) and the number of nodes, "
        "canSNPs and reference genomes.",
        epilog=examples("cansnper2 db info -d francisella_tularensis.db"),
        formatter_class=HelpFormatter,
    )
    info.add_argument("-d", "--database", metavar="DB", help="the database; required")
    add_logging_options(
        info,
        logs_dest="logs",
        logs_default=None,
        logs_help="folder for a log file; without it, db info writes no log file",
    )
    parsers["db info"] = info

    export = db_sub.add_parser(
        "export",
        help="export a database to text files",
        description="Write a database as text files: the scheme files that db "
        "import reads, or the tree in Newick format.",
        epilog=examples(
            "cansnper2 db export -d scheme.db -o export/",
            "cansnper2 db export -d scheme.db --format newick > scheme.nwk",
        ),
        formatter_class=HelpFormatter,
    )
    export.add_argument("-d", "--database", metavar="DB", help="the database; required")
    export.add_argument(
        "--format",
        dest="export_format",
        choices=database_command.supported_output,
        help="source: the tree, canSNP, reference genome and (when the database "
        "has them) genomes files that db import reads; newick: the tree only, "
        "to <outdir>/<database>_tree.nwk, or to standard output without "
        "-o; tab: former name of source (default: source)",
    )
    export.add_argument(
        "-o", "--outdir", metavar="DIR", help="folder for the exported files (default: the current folder)"
    )
    add_logging_options(export, logs_dest="logs", logs_default="logs")
    parsers["db export"] = export

    convert = db_sub.add_parser(
        "convert",
        help="convert a database or canSNP file from an earlier format",
        description="Convert a database or canSNP file written by CanSNPer2 before "
        "2.4.0 to the current format. Files from before 2.3.0 place canSNPs on "
        "concatenated genomes; their positions are moved to the reference "
        "sequences, which needs the reference files (--ref-dir).",
        epilog=examples(
            "cansnper2 db convert -d old.db --ref-dir references/ -o converted/",
            "cansnper2 db convert --cansnps old_snp.txt --ref-dir references/ \\",
            "    -o converted/",
        ),
        formatter_class=HelpFormatter,
    )
    convert.add_argument(
        "-d", "--database", metavar="DB", help="the database to convert (or give --cansnps instead)"
    )
    convert.add_argument(
        "--cansnps",
        dest="annotation",
        metavar="FILE",
        help="the canSNP source file to convert (or give --database instead)",
    )
    convert.add_argument(
        "--ref-dir",
        dest="refdir",
        metavar="DIR",
        help="folder with the reference files <genome>.fna, needed for formats "
        "before 2.3.0; missing files are downloaded",
    )
    convert.add_argument("-o", "--outdir", metavar="DIR", help="folder for the converted file; required")
    add_logging_options(convert, logs_dest="logs", logs_default="logs")
    parsers["db convert"] = convert

    modify = db_sub.add_parser(
        "modify",
        help="change the tree, canSNPs or reference genomes of a database",
        description="Change an existing database: add or replace part of the tree "
        "(--mod-file), or add canSNPs (--cansnps) or reference genomes "
        "(--references). canSNPs and reference genomes that the database "
        "already has are refused, not changed.",
        epilog=examples(
            "cansnper2 db modify -d scheme.db --mod-file new_branch.txt --parent B.Br.001",
            "cansnper2 db modify -d scheme.db --mod-file subtree.txt \\",
            "    --parent B.Br.001 --replace",
            "cansnper2 db modify -d scheme.db --cansnps added_canSNPs.txt",
        ),
        formatter_class=HelpFormatter,
    )
    modify.add_argument("-d", "--database", metavar="DB", help="the database; required")
    modify.add_argument(
        "--mod-file", dest="mod_file", metavar="FILE", help="file with the tree modifications"
    )
    modify.add_argument(
        "--parent", metavar="NODE", help="node at which the tree is changed; required with --mod-file"
    )
    modify.add_argument(
        "--replace",
        action="store_true",
        help="replace the subtree below --parent with the tree of --mod-file",
    )
    modify.add_argument(
        "--cansnps", dest="annotation", metavar="FILE", help="canSNP file with canSNPs to add"
    )
    modify.add_argument("--references", metavar="FILE", help="reference genomes file with genomes to add")
    add_logging_options(modify, logs_dest="logs", logs_default="logs")
    parsers["db modify"] = modify

    remove = db_sub.add_parser(
        "remove",
        help="remove the branch rooted at a node",
        description="Remove a node and every node below it from a database.",
        epilog=examples("cansnper2 db remove -d scheme.db --parent B.Br.001"),
        formatter_class=HelpFormatter,
    )
    remove.add_argument("-d", "--database", metavar="DB", help="the database; required")
    remove.add_argument("--parent", metavar="NODE", help="root of the branch to remove; required")
    add_logging_options(remove, logs_dest="logs", logs_default="logs")
    parsers["db remove"] = remove

    """the internal options no verb shows get their old defaults here; the options
		the verbs do show already carry their default on the option itself"""
    for name, sub in parsers.items():
        database_defaults(sub, **DB_MODES[name])
    return parsers


DB_MODES = {
    "db create": dict(create=True),
    "db import": dict(do_import=True),
    "db info": dict(info=True),
    "db export": dict(export=True),
    "db convert": dict(convert=True),
    "db modify": {},
    "db remove": dict(remove=True),
}


def add_scheme_subparsers(subparsers):
    scheme = subparsers.add_parser(
        "scheme", help="list: the published schemes; get: build the database of one"
    )
    scheme_sub = scheme.add_subparsers(dest="scheme_action", metavar="<action>", required=True)

    listing = scheme_sub.add_parser(
        "list",
        help="list the schemes and their status",
        description="List the schemes of CanSNPer2-data (or of another repository "
        "or local directory) with their status. This verb writes no log file.",
        epilog=examples("cansnper2 scheme list", "cansnper2 scheme list --branch dev"),
        formatter_class=HelpFormatter,
    )
    add_scheme_source_options(listing)
    download_defaults(listing, list_schemes=True)
    parsers = {"scheme list": listing}

    get = scheme_sub.add_parser(
        "get",
        help="fetch a scheme, build its database and download its references",
        description="Fetch the files of a scheme, check them, build the database "
        "and download its reference genomes from NCBI GenBank into ./references. "
        "Use ref get to download to another folder or from RefSeq.",
        epilog=examples(
            "cansnper2 scheme get francisella_tularensis", "cansnper2 scheme get brucella_genus --branch dev"
        ),
        formatter_class=HelpFormatter,
    )
    get.add_argument("scheme", metavar="NAME", help="name of the scheme (see scheme list)")
    get.add_argument("-d", "--database", metavar="DB", help="the database to build (default: <NAME>.db)")
    add_scheme_source_options(get)
    get.add_argument("--overwrite", action="store_true", help="replace an existing database")
    get.add_argument(
        "--no-references",
        dest="no_references",
        action="store_true",
        help="build the database without downloading its references",
    )
    get.add_argument(
        "--allow-unverified",
        dest="allow_unverified",
        action="store_true",
        help="build a verified scheme although its files cannot be checked "
        "against the checksums of its manifest",
    )
    add_local_references_option(get, "references")
    add_logging_options(get, logs_dest="logs", logs_default="logs")
    download_defaults(get)
    parsers["scheme get"] = get
    return parsers


def add_ref_subparsers(subparsers):
    ref = subparsers.add_parser("ref", help="get: download the reference genomes of a database from NCBI")
    ref_sub = ref.add_subparsers(dest="ref_action", metavar="<action>", required=True)
    get = ref_sub.add_parser(
        "get",
        help="download the reference genomes of a database",
        description="Download the reference genomes a database names from NCBI "
        "into a folder, as <genome>.fna.",
        epilog=examples(
            "cansnper2 ref get -d francisella_tularensis.db",
            "cansnper2 ref get -d scheme.db -s refseq -o references/",
        ),
        formatter_class=HelpFormatter,
    )
    get.add_argument(
        "-d", "--database", metavar="DB", help="the database whose reference genomes are downloaded; required"
    )
    get.add_argument(
        "-s", "--source", choices=["genbank", "refseq"], help="NCBI assembly database to download from"
    )
    get.add_argument("-o", "--outdir", metavar="DIR", help="folder for the reference genomes")
    add_local_references_option(get, "<outdir>")
    add_logging_options(get, logs_dest="logs", logs_default="logs")
    download_defaults(get)
    return {"ref get": get}


def add_local_references_option(parser, folder):
    parser.add_argument(
        "--local-references",
        dest="local_references",
        action="store_true",
        help="also download the representative genomes of the scheme into "
        "{folder}/local, for cansnper2 type --local-references".format(folder=folder),
    )


def download_defaults(parser, **mode):
    """every internal option of the download command with its default, plus the
    one verb of it that this subcommand runs"""
    settings = dict(
        database=None,
        scheme=None,
        branch=None,
        repo=None,
        data_dir=None,
        list_schemes=False,
        overwrite=False,
        no_references=False,
        allow_unverified=False,
        source="genbank",
        outdir="references",
        verbose=None,
        debug=None,
        local_references=False,
    )
    settings.update(mode)
    parser.set_defaults(**settings)


def build_parsers():
    """the top-level parser and a dict of every leaf subparser, keyed by its verb
    names, so a verb gets its own usage and error messages. The main commands
    are the top of the help and every verb names the subcommands or options it
    uses, which are shown in detail only in the help of that verb"""
    parser = argparse.ArgumentParser(
        prog="cansnper2",
        description="CanSNPer2 assigns bacterial genomes to lineages of a canonical SNP\n"
        "(canSNP) typing scheme, from assemblies or sequencing reads. It also\n"
        "builds and maintains scheme databases and fetches published schemes\n"
        "and their reference genomes.",
        epilog="Run cansnper2 <command> --help for the options of a command.\n"
        "Before 2.6.0 these jobs were done by the commands CanSNPer2,\n"
        "CanSNPer2-database, CanSNPer2-download and CanSNPer2-test, which are\n"
        "no longer installed; see docs/upgrading.md.\n"
        "Documentation: https://github.com/FOI-Bioinformatics/CanSNPer2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version", action="version", version="CanSNPer2 - version {version}".format(version=__version__)
    )
    subparsers = parser.add_subparsers(dest="verb", metavar="<command>", required=True)

    parsers = {"type": add_type_subparser(subparsers)}

    selftest_parser = subparsers.add_parser(
        "selftest",
        help="check the installation: Python modules and external programs",
        description="Check that the Python modules can be imported and report "
        "which aligners and read typing programs are found in PATH.",
    )
    selftest_parser.add_argument(
        "--debug", action="store_true", help="print the traceback of an unexpected error"
    )

    parsers.update(add_db_subparsers(subparsers))
    parsers.update(add_scheme_subparsers(subparsers))
    parsers.update(add_ref_subparsers(subparsers))
    return parser, parsers


def run_command(args, parsers):
    """call the same run function the matching old command calls"""
    if args.verb == "type":
        typing_command.run_command(args, parsers["type"])
    elif args.verb == "selftest":
        selftest_command.run_command()
    elif args.verb == "db":
        database_command.run_command(args)
    elif args.verb == "scheme":
        key = "scheme " + args.scheme_action
        download_command.run_command(args, parsers[key])
    elif args.verb == "ref":
        download_command.run_command(args, parsers["ref get"])


def command():
    """the entry point of cansnper2: one parser tree, one verb dispatch"""
    parser, parsers = build_parsers()
    args = parser.parse_args()
    run_command(args, parsers)


def main():
    """the entry point of the command: command() with the handler of unexpected errors"""
    return run_main(command)


if __name__ == "__main__":
    main()
