import os

from sv_util.settings import *

def sniffles(bam_file, vcf_file, reference_path):
    # threads: -t
    # Minimum number of reads that support a SV: -s
    os.system(sniffles_path + " -t 32 -s 2 -m " + bam_file + " -v " + vcf_file
                + " >/dev/null 2>&1")


def delly(bam_file, vcf_file, reference_path):
    os.system(delly_path + " call -g " + reference_path + " " + bam_file
                + " -o " + vcf_file + ".bcf >/dev/null 2>&1")

    if os.path.exists(vcf_file + ".bcf"):
        os.system(bcf_tools_path + " view " + vcf_file + ".bcf > " + vcf_file)

def gridss(bam_file, vcf_file, reference_path):
    curr = os.getcwd()
    if os.path.exists(bam_file + ".gridss.work"):
        os.system('rm -rf ' + bam_file + ".gridss.work")
    os.mkdir(bam_file + ".gridss.work")
    os.chdir(bam_file + ".gridss.work")
    os.system(gridss_path + " -t 32 -r " + reference_path + " -o " + vcf_file + " -a assembly.bam " + bam_file + 
              " >/dev/null 2>&1")
    os.chdir(curr)

def manta(bam_file, vcf_file, reference_path):
    # prepare manta
    if os.path.exists(vcf_file + ".manta"):
        os.system('rm -rf ' + vcf_file + ".manta")
    os.mkdir(vcf_file + ".manta")
    #manta_path = "~/miniconda3/envs/manta/share/manta-1.6.0-1/bin/configManta.py"
    manta_path = "~/workspace/manta/install/bin/configManta.py >/dev/null 2>&1"
    os.system("python2 "+ manta_path + " --referenceFasta " + reference_path + " --bam " + bam_file + " --runDir " + vcf_file + ".manta" )
    # actually run manta
    os.system("python2 " + vcf_file + ".manta/runWorkflow.py -j 32 -m local >" + vcf_file + ".manta/workflow.stdout.log.txt 2>" + vcf_file + ".manta/workflow.stderr.log.txt" )

    os.system("cp " + vcf_file + ".manta/results/variants/diploidSV.vcf.gz" + " " + vcf_file + ".gz")
    os.system("gunzip -f " + vcf_file + ".gz")
