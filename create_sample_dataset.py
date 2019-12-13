from MA import *
import os
import textwrap
import time
import json
import compare_callers
import random
import subprocess
import random

supporting_nt = 10**6
coverage = 1

global_prefix = "C:/MAdata/"

# AKFIX
"""Markus @ Zeus""" 
svdb_dir = "/MAdata/sv_datasets/" # AKFIX
survivor = "~/workspace/SURVIVOR/Debug/SURVIVOR simreads "  
genome_dir = "/MAdata/genome/human"
survivor_error_profile_dir = "~/workspace/SURVIVOR/"
OS_is_MSWIN = False

"""Arne @ home """
# survivor = global_prefix + "tools/Survivor.exe simreads " # Arne @ desktop at home
# svdb_dir = global_prefix + "sv_datasets/" 
# genome_dir = global_prefix + "genome/GRCh38.p12-chr1"
# survivor_error_profile_dir = global_prefix + "tools/"
# OS_is_MSWIN = True


def create_illumina_reads_dwgsim(sequenced_genome_pack, ref_pack, sequenced_genome_path, database, reads_folder,
                                 json_info_file, coverage, name, read_length):
    json_info_file["read_length"] = read_length
    reads1 = reads_folder + name + ".bwa.read1.fastq.gz"
    reads2 = reads_folder + name + ".bwa.read2.fastq.gz"
    json_info_file["fasta_file"] = reads1
    json_info_file["fasta_file_mate"] = reads2

    print("\tdwgsim...")
    dwgsim_instances = []
    file_names1 = ""
    file_names2 = ""
    # we actually need to seed dwgsim otherwise each instance will create the same data (it uses the current time)
    seed = random.randrange(0, 2**6)
    num_instances = 32
    for idx in range(num_instances):
        dwgsim = "/usr/home/markus/workspace/DWGSIM/dwgsim"
        #command = [dwgsim, "-1", str(read_length), "-2", str(read_length), "-e", "from 0.000 to 0.000 by 0.000", \
        #           "-E", "from 0.000 to 0.000 by 0.000", "-y", "0", "-C", str(coverage/32), "-r", "0", "-o", \
        #           "1", "-z", str(seed), sequenced_genome_path, reads_folder + "part_" + str(idx) + name]
        command = [dwgsim, "-1", str(read_length), "-2", str(read_length), "-C", str(coverage/num_instances), "-r", "0", "-o", \
                   "1", "-z", str(seed), sequenced_genome_path, reads_folder + "part_" + str(idx) + name]
        dwgsim_instances.append(subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        seed += 42 # increment the seed so that we have some different number every time
        file_names1 += reads_folder + "part_" + str(idx) + name + ".bwa.read1.fastq.gz "
        file_names2 += reads_folder + "part_" + str(idx) + name + ".bwa.read2.fastq.gz "
    for proc in dwgsim_instances:
        proc.wait()

    print("\tcat...")
    os.system("zcat " + file_names1 + " > " + reads1)
    os.system("zcat " + file_names2 + " > " + reads2)

    print("\tinserting into db...")
    inserter = ReadInserter(database, name, ref_pack)
    json_info_file["seq_id"] = inserter.sequencer_id
    f_path_vec_1 = libMA.filePathVector([libMA.path(reads1)])
    f_path_vec_2 = libMA.filePathVector([libMA.path(reads2)])
    inserter.insert_paired_fasta_files(ParameterSetManager(), f_path_vec_1, f_path_vec_2)
    print("\tdone")

def create_reads_survivor(sequenced_genome_pack, ref_pack, sequenced_genome_path, database, reads_folder, 
                          json_info_file, coverage, name, error_profile, technology):
    json_info_file["error_profile"] = error_profile
    json_info_file["technology"] = technology
    reads1 = reads_folder + name + ".fasta"
    json_info_file["fasta_file"] = reads1

    command = survivor + sequenced_genome_path + " " + error_profile + " " + str(coverage) + " " \
              + reads1
    print("Command:", command)

    if OS_is_MSWIN:
        os.system(command)
    else:
        os.system(command + " >/dev/null 2>&1")

    print("\tinserting into db...")
    inserter = ReadInserter(database, name, ref_pack)
    json_info_file["seq_id"] = inserter.sequencer_id
    f_path_vec = libMA.filePathVector([libMA.path(reads1)])
    inserter.insert_fasta_files(ParameterSetManager(), f_path_vec)
    print("\tdone")

def sv_deletion(sv_inserter, position, sv_size):
    sv_inserter.insert_call(SvCall(position, position + sv_size, 0, 0, False, supporting_nt, coverage))

def sv_duplication(sv_inserter, position, sv_size):
    sv_inserter.insert_call(SvCall(position + sv_size, position, 0, 0, False, supporting_nt, coverage))

def sv_inversion(sv_inserter, position, sv_size):
    sv_inserter.insert_call(SvCall(position + sv_size, position, 0, 0, True, supporting_nt, coverage))
    sv_inserter.insert_call(SvCall(position, position + sv_size, 0, 0, True, supporting_nt, coverage))

def sv_translocation(sv_inserter, position, sv_size, gap_size):
    assert gap_size < sv_size # sanity check
    assert gap_size > 10 # otherwise (1) & (2) combine & we need a special case for that
    start_a = position
    end_a = position + int((sv_size - gap_size) / 2)
    start_b = end_a + gap_size
    end_b = position + sv_size
    sv_inserter.insert_call(SvCall(start_a, start_b, 0, 0, False, supporting_nt, coverage))
    sv_inserter.insert_call(SvCall(end_b, end_a, 0, 0, False, supporting_nt, coverage)) # (1)
    sv_inserter.insert_call(SvCall(start_b - 1, start_a + 1, 0, 0, False, supporting_nt, coverage)) # (2)
    sv_inserter.insert_call(SvCall(end_a - 1, end_b + 1, 0, 0, False, supporting_nt, coverage))

def sv_insertion(sv_inserter, position, sv_size):
    call = SvCall(position, position + 1, 0, 0, False, supporting_nt, coverage)

    ins = ""
    for _ in range(sv_size):
        ins += random.choice(["a", "c", "g", "t"])

    call.inserted_sequence = NucSeq(ins)
    sv_inserter.insert_call(call)

##
# sv_func signature: def sv_func(sv_inserter, position, sv_size)
#
def separate_svs(pack, database, json_info_file, sv_func, sv_size, sv_margin, chromosome=None):
    json_info_file["sv_size"] = sv_size
    json_info_file["sv_margin"] = sv_margin
    json_info_file["sv_func"] = sv_func[0].__name__

    # -1 since there are no related sv jumps
    sv_inserter = SvCallInserter(database, json_info_file["name"] + "_simulated_sv", "the sv's that were simulated", -1)

    def x(s, l):
        for pos in range(s + sv_margin, s + l - sv_margin, sv_size + sv_margin):
            sv_func[0](sv_inserter, pos, sv_size, *sv_func[1])
    if chromosome is None:
        for s, l in zip(pack.contigStarts(), pack.contigLengths()):
            x(s, l)
    else:
        s = pack.contigStarts()[pack.id_of_sequence(chromosome)]
        l = pack.contigLengths()[pack.id_of_sequence(chromosome)]
        x(s, l)

    return sv_inserter.sv_caller_run_id

def no_svs(pack, database, json_info_file):
    sv_inserter = SvCallInserter(database, json_info_file["name"] + "_simulated_sv", "the sv's that were simulated", -1)
    return sv_inserter.sv_caller_run_id


##
# create_svs_func signature: def create_svs_func(pack, database, json_info_file) -> caller_id
# create_reads_funcs is a list of functions with the signature: 
#       def create_reads_func(sequenced_genome_pack, sequenced_genome_path, database, reads_folder,
#                             json_info_file, coverage, name)
#
def create_dataset(reference_path, # dir with reference 
                   dataset_name, create_svs_funcs,
                   create_reads_funcs, coverages, chromosome=None):
    ref_pack = Pack()
    ref_pack.load(reference_path + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(reference_path + "/ma/genome")
    if not chromosome is None:
        idx = ref_pack.id_of_sequence(chromosome)
        print("picked contig:", idx)
        print("startindex, endindex: ", ref_pack.start_of_sequence_id(idx), 
                                      ref_pack.start_of_sequence_id(idx) + ref_pack.length_of_sequence(chromosome))

    os.mkdir(svdb_dir + dataset_name) # this throws an error if the dataset already exists
    os.mkdir(svdb_dir + dataset_name + "/reads")
    os.mkdir(svdb_dir + dataset_name + "/genomes")
    os.mkdir(svdb_dir + dataset_name + "/alignments")
    os.mkdir(svdb_dir + dataset_name + "/calls")

    json_info_file = {
        "reference_path": reference_path,
        "datasets": [],
        "version": 1
    }

    print("creating SV-DB...")
    start = time.time()

    # create the sv_db
    database = SV_DB(svdb_dir + dataset_name + "/svs.db", "create")
    print(time.time() - start, "seconds")

    for create_svs_func, sv_func_name, create_svs_funcs_params in create_svs_funcs:
        seq_gen_path = svdb_dir + dataset_name + "/genomes/sequenced_genome_" + sv_func_name
        print("creating", sv_func_name, "dataset ...")
        start = time.time()
        # create the svs
        json_info_file_dataset_sub = {
            "func_name": create_svs_func.__name__,
            "create_reads_funcs": [],
            "sequenced_genome_path": seq_gen_path,
            "name": sv_func_name
        }
        caller_id = create_svs_func(ref_pack, database, json_info_file_dataset_sub, *create_svs_funcs_params)
        json_info_file_dataset_sub["ground_truth"] = caller_id
        print(time.time() - start, "seconds")

        print("creating sequenced genome...")
        start = time.time()

        # save the sequenced genome
        
        print("\treconstructing sequenced genome")
        seq_pack = database.reconstruct_sequenced_genome(ref_pack, caller_id)
        if not chromosome is None:
            seq_pack_ = Pack()
            seq_pack_.append(seq_pack.contigNames()[ref_pack.id_of_sequence(chromosome)], "no_description_given", 
                             NucSeq(seq_pack.contigSeqs()[ref_pack.id_of_sequence(chromosome)]))
            seq_pack = seq_pack_
        print("\tstoring sequenced genome")
        seq_pack.store(seq_gen_path)
        with open(seq_gen_path + ".fasta", "w") as fasta_out:
            for name, sequence in zip(seq_pack.contigNames(), seq_pack.contigSeqs()):
                print("writing:", name, "len:", len(sequence))
                fasta_out.write(">")
                fasta_out.write(name)
                fasta_out.write("\n")
                idx = 0
                while idx < len(sequence):
                    fasta_out.write(sequence[idx:idx+50])
                    fasta_out.write("\n")
                    idx += 50
        print(time.time() - start, "seconds")

        if not chromosome is None:
            print("\tvalidating created genome...")
            seeding_module = BinarySeeding(ParameterSetManager())
            q = seq_pack.extract_forward_strand_n()
            segments = seeding_module.execute(fm_index, q)
            num_nts = [0]*len(ref_pack.contigLengths())
            for seeds in segments.extract_seeds(fm_index, 1, 16, len(q), True):
                pos = ref_pack.seq_id_for_pos(seeds.start_ref)
                num_nts[pos] += seeds.size
            highest_idx = 0
            print("id\tcoverage\tnt generated")
            for idx, (num_nt, size) in enumerate(zip(num_nts, ref_pack.contigLengths())):
                if round(100*num_nt/size, 3) >= 3:
                    print(idx, str(round(100*num_nt/size, 3)) + "%", num_nt, sep='\t')
                    if num_nt > num_nts[highest_idx]:
                        highest_idx = idx
            if highest_idx != ref_pack.id_of_sequence(chromosome):
                print("WARNING picked chromosome seems not to be the one that was generated")


        print("\tdone")
        print("creating reads...")
        start = time.time()

        # create the reads
        for create_reads_func, name, create_reads_args in create_reads_funcs:
            for coverage in coverages:
                name_c = sv_func_name + "-" + name + "-" + str(coverage).zfill(2) + "x"
                print(name_c, "...")
                json_info_file_sub = {
                    "func_name": create_reads_func.__name__,
                    "name": name_c,
                    "coverage": coverage
                }
                create_reads_func(seq_pack, ref_pack, seq_gen_path + ".fasta",
                                database, svdb_dir + dataset_name + "/reads/", json_info_file_sub,
                                coverage, name_c, *create_reads_args)
                json_info_file_dataset_sub["create_reads_funcs"].append(json_info_file_sub)
        json_info_file["datasets"].append(json_info_file_dataset_sub)
        print(time.time() - start, "seconds")

    # save the info.json file
    with open(svdb_dir + dataset_name + "/info.json", "w") as json_out:
        json.dump(json_info_file, json_out)

    print("done creating dataset:")
    print(json_info_file)


if __name__ == "__main__":
    survivor_error_profile_pac_b = survivor_error_profile_dir + "HG002_Pac_error_profile_bwa.txt"
    survivor_error_profile_ont = survivor_error_profile_dir + "NA12878_nano_error_profile_bwa.txt"

    create_dataset(genome_dir + "/GRCh38.p12-chr1",
                   "minimal",
                   [
                    ( separate_svs, "del-0100", ( (sv_deletion, tuple()), 100, 1000 ) ),
                    ],
                   [(create_illumina_reads_dwgsim, "ill_250", (250,)),],
                   [10])

    #create_dataset(genome_dir + "/GRCh38.p12-chr1-large",
    #               "comprehensive",
    #               [( separate_svs, "del-0100", ( (sv_deletion, tuple()), 100, 500 ) ),
    #                ( separate_svs, "inv-0100", ( (sv_inversion, tuple()), 100, 500 ) ),
    #                ( separate_svs, "dup-0100", ( (sv_duplication, tuple()), 100, 500 ) ),
    #                ( separate_svs, "trans-0100", ( (sv_translocation, (25,)), 100, 500 ) ),
    #                ( separate_svs, "ins-0100", ( (sv_insertion, tuple()), 100, 500 ) ),

    #                ( separate_svs, "del-0250", ( (sv_deletion, tuple()), 250, 1000 ) ),
    #                ( separate_svs, "inv-0250", ( (sv_inversion, tuple()), 250, 1000 ) ),
    #                ( separate_svs, "dup-0250", ( (sv_duplication, tuple()), 250, 1000 ) ),
    #                ( separate_svs, "trans-0250", ( (sv_translocation, (50,)), 250, 1000 ) ),
    #                ( separate_svs, "ins-0250", ( (sv_insertion, tuple()), 250, 1000 ) ),

    #                ( separate_svs, "del-1000", ( (sv_deletion, tuple()), 1000, 5000 ) ),
    #                ( separate_svs, "inv-1000", ( (sv_inversion, tuple()), 1000, 5000 ) ),
    #                ( separate_svs, "dup-1000", ( (sv_duplication, tuple()), 1000, 5000 ) ),
    #                ( separate_svs, "trans-1000", ( (sv_translocation, (200,)), 1000, 5000 ) ),
    #                ( separate_svs, "ins-1000", ( (sv_insertion, tuple()), 1000, 5000 ) )],
    #               [(create_illumina_reads_dwgsim, "ill_250", (250,)),
    #                (create_illumina_reads_dwgsim, "ill_150", (150,)),
    #                (create_illumina_reads_dwgsim, "ill_100", (100,)),
    #                (create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb")),
    #                (create_reads_survivor, "ont", (survivor_error_profile_ont, "ont"))
    #                ],
    #               [5, 10, 25, 50])

    #chrom = "CM000679.2" # Chromosome 17
    #create_dataset("/MAdata/genome/human/GRCh38.p12",
    #               "del_human",
    #               [( separate_svs, "del-100", ( (sv_deletion, tuple()), 100, 500 ) ),],
    #               [(create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb"))],
    #               [25])

    #for prefix, func in [ ("del", sv_deletion), ("ins", sv_insertion), ("dup",sv_duplication), ("inv", sv_inversion) ]:
    #    chrom = "CM000663.2" # Chromosome 1
    #    create_dataset("/MAdata/genome/human/GRCh38.p12",
    #                prefix + "_human",
    #                [( separate_svs, prefix + "-0250", ( (func, tuple()), 250, 1000, chrom ) ),
    #                    ( separate_svs, prefix + "-1000", ( (func, tuple()), 1000, 5000, chrom ) ),],
    #                [(create_illumina_reads_dwgsim, "ill_250", (250,)),
    #                    (create_illumina_reads_dwgsim, "ill_150", (150,)),
    #                    (create_illumina_reads_dwgsim, "ill_100", (100,)),
    #                    (create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb"))],
    #                [5, 10, 25],
    #                chrom)
    #
    #chrom = "CM000663.2" # Chromosome 1
    #create_dataset("/MAdata/genome/human/GRCh38.p12",
    #               "tra_human",
    #               [( separate_svs, "tra-10000", ( (sv_translocation, (200,)), 1000, 5000, chrom ) ),],
    #               [(create_illumina_reads_dwgsim, "ill_250", (250,)),
    #                (create_illumina_reads_dwgsim, "ill_150", (150,)),
    #                (create_illumina_reads_dwgsim, "ill_100", (100,)),
    #                (create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb"))],
    #               [5, 10, 25],
    #               chrom)

    #chrom = "CM000683.2" # Chromosome 21 -> shortest chromosome
    #create_dataset("/MAdata/genome/human/GRCh38.p12",
    #               "comprehensive_human",
    #               [( separate_svs, "del-0250", ( (sv_deletion, tuple()), 250, 1000, chrom ) ),
    #                ( separate_svs, "inv-0250", ( (sv_inversion, tuple()), 250, 1000, chrom ) ),
    #                ( separate_svs, "dup-0250", ( (sv_duplication, tuple()), 250, 1000, chrom ) ),
    #                ( separate_svs, "trans-0250", ( (sv_translocation, (50,)), 250, 1000, chrom ) ),
    #                ( separate_svs, "ins-0250", ( (sv_insertion, tuple()), 250, 1000, chrom ) ),
    #                ( separate_svs, "del-1000", ( (sv_deletion, tuple()), 1000, 5000, chrom ) ),
    #                ( separate_svs, "inv-1000", ( (sv_inversion, tuple()), 1000, 5000, chrom ) ),
    #                ( separate_svs, "dup-1000", ( (sv_duplication, tuple()), 1000, 5000, chrom ) ),
    #                ( separate_svs, "trans-1000", ( (sv_translocation, (200,)), 1000, 5000, chrom ) ),
    #                ( separate_svs, "ins-1000", ( (sv_insertion, tuple()), 1000, 5000, chrom ) )],
    #               [(create_illumina_reads_dwgsim, "ill_250", (250,)),
    #                (create_illumina_reads_dwgsim, "ill_150", (150,)),
    #                (create_illumina_reads_dwgsim, "ill_100", (100,)),
    #                (create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb"))],
    #               [5, 10, 25],
    #               chrom)

    """
    for sv_size in [100, 150, 500]:
        for sv_type in [sv_deletion, sv_duplication, sv_inversion, sv_insertion]:
            dataset_name = "separate_random_" + sv_type.__name__ + "_" + str(sv_size) + "nt"
            create_dataset("/MAdata/genome/random",
                        dataset_name,
                        ( create_separate_svs, ( (sv_type, tuple()), sv_size, 1000 ) ),
                        [(create_illumina_reads_dwgsim, "illumina-150nt", (150,)),
                            (create_illumina_reads_dwgsim, "illumina-250nt", (250,)),
                            #(create_reads_survivor, "ont", (survivor_error_profile_ont, "ont")),
                            (create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb"))],
                        [5, 10, 25, 50])
            analyze_sample_dataset(dataset_name)
    """