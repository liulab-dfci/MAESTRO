import os,sys
import snakemake.report

import argparse as ap

from MAESTRO.scATAC_utility import is_gzip, get_fastqlist, ENV_PATH, SCRIPT_PATH, RSCRIPT_PATH

def CommandLineParser():
    parser = ap.ArgumentParser(description = "Generate scATAC result report. ")

    group_input = parser.add_argument_group("Input arguments")
    group_input.add_argument("--platform", dest = "platform", default = "10x-genomics",
        choices = ["10x-genomics", "sci-ATAC-seq", "microfluidic"],
        help = "Platform of single cell ATAC-seq. DEFAULT: 10x-genomics.")
    group_input.add_argument("--input-format", dest = "input_format", default = "fastq",
        choices = ["fastq", "bam", "fragments"],
        help = "The format of input files. DEFAULT: fastq.")
    group_input.add_argument("--input-path", dest = "input_path", type = str, default = "",
        help = "Directory where input files are stored")
    group_input.add_argument("--species", dest = "species", default = "GRCh38",
        choices = ["GRCh38", "GRCm38"], type = str,
        help = "Species (GRCh38 for human and GRCm38 for mouse). DEFAULT: GRCh38.")

    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-d", "--directory", dest = "directory", default = "MAESTRO",
        help = "Path to the directory where the result file shall be stored. DEFAULT: MAESTRO.")
    group_output.add_argument("--outprefix", dest = "outprefix", default = "10x-genomics",
        help = "Prefix of output files. DEFAULT: MAESTRO.")

    return parser.parse_args()


def main():

    #SCRIPT_PATH = os.path.dirname(__file__)

    myparser = CommandLineParser()
    directory = myparser.directory
    outpre = myparser.outprefix
    input_path = myparser.input_path
    species = myparser.species
    platform = myparser.platform
    input_format = myparser.input_format

    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass

    batch_peakcluster_link = snakemake.report.data_uri_from_file("Result/Analysis/Batch/all_samples_cluster.png")
    batch_cluster_regulator_file = "Result/Analysis/Batch/all_samples.PredictedTFTop10.txt"

    if os.path.exists("Result/Analysis/Batch/all_samples_annotated.png"):
        batch_rpannodisplay = "inline"
        batch_rpannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/Batch/all_samples_annotated.png")
    else:
        batch_rpannodisplay = "none"
        batch_rpannotate_link = ""
    if os.path.exists("Result/Analysis/Batch/all_samples_CistromeTop_annotated.png"):
        batch_caannodisplay = "inline"
        batch_caannotate_link = snakemake.report.data_uri_from_file("Result/Analysis/Batch/all_samples_CistromeTop_annotated.png")
    else:
        batch_caannodisplay = "none"
        batch_caannotate_link = ""
    if os.path.exists("Result/Analysis/Batch/all_samples_MS4A1_genetrack.png"):
        ms4a1display = "inline"
        ms4a1track_link = snakemake.report.data_uri_from_file("Result/Analysis/Batch/all_samples_MS4A1_genetrack.png")
    else:
        ms4a1display = "none"
        ms4a1track_link = ""

    if os.path.exists("Result/Analysis/Batch/all_samples_CD3D_genetrack.png"):
        cd3ddisplay = "inline"
        cd3dtrack_link = snakemake.report.data_uri_from_file("Result/Analysis/Batch/all_samples_CD3D_genetrack.png")
    else:
        cd3ddisplay = "none"
        cd3dtrack_link = ""

    batch_td_list = []
    for line in open(batch_cluster_regulator_file,"r").readlines():
        if not line.startswith("Cluster"):
            batch_items = line.strip().split("\t")
            batch_items_str_list = ["                                                            <td>" + i + "</td>" for i in batch_items]
            batch_items_str = "                                                        <tr>\n" + "\n".join(batch_items_str_list) + "\n                                                        </tr>"
            batch_td_list.append(batch_items_str)
    batch_td_str = "\n".join(batch_td_list)

    if input_format != "fragments":
        if input_format == "fastq":
            inputformat = "FASTQ Path"
        else:
            inputformat = "BAM Path"
        report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_template.html")
        report_html_temp = open(report_html_tempfile, "r").read()

        report_html_instance = report_html_temp % {"outprefix":outpre, "inputformat":inputformat, "inputpath":input_path, "species":species,
        "platform":platform, "peakcluster": batch_peakcluster_link, "rpannodisplay": batch_rpannodisplay, "rpannotate":batch_rpannotate_link, "regtable":batch_td_str,
        "caannodisplay": batch_caannodisplay, "caannotate": batch_caannotate_link, 'cd3ddisplay': cd3ddisplay, 'ms4a1display': ms4a1display,
        'cd3dtrack': cd3dtrack_link, "ms4a1track": ms4a1track_link}
    else:
        inputformat = "Fragment Path"
        report_html_tempfile = os.path.join(SCRIPT_PATH, "html", "scATAC_nomappability_template.html")
        report_html_temp = open(report_html_tempfile, "r").read()

        report_html_instance = report_html_temp % {"outprefix":outpre, "inputformat":inputformat, "inputpath":input_path, "species":species,
        "platform":platform, "peakcluster":batch_peakcluster_link, "rpannodisplay": batch_rpannodisplay, "rpannotate":batch_rpannotate_link, "regtable":batch_td_str,
        "caannodisplay": batch_caannodisplay, "caannotate": batch_caannotate_link, 'cd3ddisplay': cd3ddisplay, 'ms4a1display': ms4a1display,
        'cd3dtrack': cd3dtrack_link, "ms4a1track": ms4a1track_link}

    report_html_instancefile = os.path.join(directory, "all_samples_scATAC_report.html")
    outf = open(report_html_instancefile,"w")
    outf.write(report_html_instance)
    outf.close()


if __name__ == "__main__":
    main()
