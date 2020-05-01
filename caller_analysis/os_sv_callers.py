def sniffles(bam_file, vcf_file):
    # threads: -t
    # Minimum number of reads that support a SV: -s
    os.system("~/workspace/Sniffles/bin/sniffles-core-1.0.8/sniffles -t 32 -s 2 -m " + bam_file + " -v " + vcf_file
                + " >/dev/null 2>&1")

def pbHoney(bam_file, vcf_file):
    os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py pie -n 32 -o " + vcf_file + ".pie.bam " + bam_file +
                " " + json_dict["reference_path"] + "/blasr/genome.fasta")

    os.system("~/workspace/samtools/samtools sort -@ 32 -m 1G " + vcf_file + ".pie.bam > " + 
                vcf_file + ".pie.sorted.bam")
    os.system("~/workspace/samtools/samtools index " + vcf_file + ".pie.sorted.bam > " + 
                vcf_file + ".pie.sorted.bam.bai")

    with open(vcf_file + ".hon.tails", "w") as out_file:
        out_file.write("##fileformat=VCFv4.0\n")
        out_file.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
    os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py tails -o " + vcf_file + ".hon.tails " + 
                vcf_file + ".pie.sorted.bam")

    with open(vcf_file + ".hon.spots.spots", "w") as out_file:
        out_file.write("##fileformat=VCFv4.0\n")
        out_file.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
    #os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py spots -n 32 -o " + vcf_file + 
    #          ".hon.spots --reference " + json_dict["reference_path"] + "/blasr/genome.fasta " + bam_file )
    os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py spots --consensus None -n 32 -o " + vcf_file + 
                ".hon.spots " + bam_file )

    with open(vcf_file, "w") as out_file:
        with open(vcf_file + ".hon.spots.spots", "r") as in_file:
            for line in in_file:
                out_file.write(line)
        with open(vcf_file + ".hon.tails", "r") as in_file:
            for line in in_file:
                if line[0] != "#":
                    out_file.write(line)


def delly(bam_file, vcf_file):
    os.system("~/workspace/delly/delly call -g " + json_dict["reference_path"] + "/fasta/genome.fna " + bam_file
                + " -o " + vcf_file + ".bcf ") #>/dev/null 2>&1

    os.system("~/workspace/bcftools/bcftools-1.9/bcftools view " + vcf_file + ".bcf > " + vcf_file)

def manta(bam_file, vcf_file):
    # prepare manta
    os.system("rm -r " + vcf_file + ".manta") # clean up folder
    os.system("python2 ~/workspace/manta/manta-1.5.0.centos6_x86_64/bin/configManta.py --referenceFasta " + json_dict["reference_path"] + "/fasta/genome.fna --bam " + bam_file + " --runDir " + vcf_file + ".manta" )
    # actually run manta
    os.system("python2 " + vcf_file + ".manta/runWorkflow.py -j 32 -m local" )

    os.system("cp " + vcf_file + ".manta/results/variants/diploidSV.vcf.gz" + " " + vcf_file + ".gz")
    os.system("gunzip -f " + vcf_file + ".gz")

def smoove(bam_file, vcf_file):
    docker = True
    if docker:
        bam_folder = bam_file[:bam_file.rfind("/")]
        bam_filename = bam_file[bam_file.rfind("/")+1:]
        vcf_folder = vcf_file[:vcf_file.rfind("/")]
        print(vcf_folder)
        os.system( "rm -f " + vcf_folder + "/*.disc.*" )
        vcf_filename = vcf_file[vcf_file.rfind("/")+1:]
        # set: -e SMOOVE_KEEP_ALL=KEEP otherwise homozygous variants are discarded
        s = "docker run -v " + bam_folder + ":/bam_folder/ -v " + json_dict["reference_path"] + \
            "/fasta:/genome_folder/ -v " + vcf_folder + \
            ":/vcf_folder/ -e SMOOVE_KEEP_ALL=KEEP -it brentp/smoove smoove call -o /vcf_folder/ --name " + \
            vcf_filename + \
            " --excludechroms dummy --fasta /genome_folder/genome.fna -p 32 -S 2 --genotype /bam_folder/" + \
            bam_filename #+ " > /dev/null"
        #print(s)
        os.system( s )
    else:
        os.system( "smoove call -o " + vcf_file + " --noextrafilters --fasta " + json_dict["reference_path"] 
                    + "/fasta/genome.fna -p 32 --genotype " + bam_file )
    os.system("gunzip -f " + vcf_file + "-smoove.genotyped.vcf.gz")
    os.system("mv " + vcf_file + "-smoove.genotyped.vcf " + vcf_file)

def pbSv(bam_file, vcf_file):
    os.system("~/miniconda2/bin/pbsv discover " + bam_file + " " + vcf_file + ".svsig.gz")
    os.system("~/miniconda2/bin/pbsv call " + json_dict["reference_path"] + "/fasta/genome.fa " + 
                vcf_file + ".svsig.gz " + vcf_file)