from MA import *
import os
import textwrap
import time
import json
import compare_callers

def create_illumina_reads_dwgsim(sequenced_genome_pack, sequenced_genome_path, database, reads_folder, json_info_file,
                                 coverage, name, read_length):
    json_info_file["read_length"] = read_length
    reads1 = reads_folder + name + ".bwa.read1.fastq.gz"
    reads2 = reads_folder + name + ".bwa.read2.fastq.gz"
    json_info_file["fasta_file"] = reads1
    json_info_file["fasta_file_mate"] = reads2

    dwgsim = "~/workspace/DWGSIM/dwgsim"
    command = dwgsim + " -1 " + str(read_length) + " -2 " + str(read_length) + \
        " -C " + str(coverage) + " " + sequenced_genome_path + " " + reads_folder + name
    os.system(command + " >/dev/null 2>&1")
    
    reader = PairedFileReader(ParameterSetManager(), 
                              libMA.filePathVector([libMA.path(reads1)]),
                              libMA.filePathVector([libMA.path(reads2)]))

    counter = 0
    inserter = ReadInserter(database, name)
    json_info_file["seq_id"] = inserter.sequencer_id
    while not reader.is_finished():
        reads = reader.execute()
        reads[0].name += "_paired_read_prim_" + str(counter)
        reads[1].name += "_paired_read_mate_" + str(counter)
        inserter.insert_paired_read(reads[0], reads[1])
        counter += 1

def create_reads_survivor(sequenced_genome_pack, sequenced_genome_path, database, reads_folder, json_info_file,
                          coverage, name, error_profile, technology):
    json_info_file["error_profile"] = error_profile
    json_info_file["technology"] = technology
    reads1 = reads_folder + name + ".fasta"
    json_info_file["fasta_file"] = reads1

    survivor = "~/workspace/SURVIVOR/Debug/SURVIVOR simreads "
    command = survivor + sequenced_genome_path + " " + error_profile + " " + str(coverage) + " " \
              + reads1

    os.system(command + " >/dev/null 2>&1")
    
    reader = FileReader(ParameterSetManager(), libMA.path(reads1))

    counter = 0
    inserter = ReadInserter(database, name)
    json_info_file["seq_id"] = inserter.sequencer_id
    while not reader.is_finished():
        read = reader.execute()
        read.name += "_read_" + str(counter)
        inserter.insert_read(read)
        counter += 1

def sv_deletion(sv_inserter, position, sv_size):
    sv_inserter.insert_call(SvCall(position, position + sv_size, 1, 1, False, float('inf')))

def sv_duplication(sv_inserter, position, sv_size):
    sv_inserter.insert_call(SvCall(position + sv_size, position, 1, 1, False, float('inf')))

def sv_inversion(sv_inserter, position, sv_size):
    sv_inserter.insert_call(SvCall(position + sv_size, position, 1, 1, True, float('inf')))
    sv_inserter.insert_call(SvCall(position, position + sv_size, 1, 1, True, float('inf')))

def sv_translocation(sv_inserter, position, sv_size, gap_size):
    assert gap_size < sv_size # sanity check
    assert gap_size > 30 # otherwise (1) & (2) combine & we need a special case for that
    start_a = position
    end_a = position + (sv_size - gap_size) / 2
    start_b = end_a + gap_size
    end_b = position + sv_size
    sv_inserter.insert_call(SvCall(start_a, start_b, 1, 1, False, float('inf')))
    sv_inserter.insert_call(SvCall(end_b, end_a, 1, 1, False, float('inf'))) # (1)
    sv_inserter.insert_call(SvCall(start_b - 1, start_a + 1, 1, 1, False, float('inf'))) # (2)
    sv_inserter.insert_call(SvCall(end_a - 1, end_b + 1, 1, 1, False, float('inf')))

def sv_insertion(sv_inserter, position, sv_size):
    call = SvCall(position, position + 1, 1, 1, False, float('inf'))

    ins = ""
    for _ in range(sv_size):
        ins += random.choice(["a", "c", "g", "t"])

    call.inserted_sequence = NucSeq(ins)
    sv_inserter.insert_call(call)

##
# sv_func signature: def sv_func(sv_inserter, position, sv_size)
#
def create_separate_svs(pack, database, json_info_file, sv_func, sv_size, sv_margin):
    json_info_file["sv_size"] = sv_size
    json_info_file["sv_margin"] = sv_margin
    json_info_file["sv_func"] = sv_func[0].__name__

    # -1 since there are no related sv jumps
    sv_inserter = SvCallInserter(database, "simulated sv", "the sv's that were simulated", -1)

    for s, l in zip(pack.contigStarts(), pack.contigLengths()):
        for pos in range(s + sv_margin, s + l - sv_margin, sv_size + sv_margin):
            sv_func[0](sv_inserter, pos, sv_size, *sv_func[1])
    return sv_inserter.sv_caller_run_id
def no_svs(pack, database, json_info_file):
    return 0


##
# create_svs_func signature: def create_svs_func(pack, database, json_info_file) -> caller_id
# create_reads_funcs is a list of functions with the signature: 
#       def create_reads_func(sequenced_genome_pack, sequenced_genome_path, database, reads_folder,
#                             json_info_file, coverage, name)
#
def create_dataset(reference_path, dataset_name, create_svs_func,
                   create_reads_funcs, coverages):
    os.mkdir("/MAdata/sv_datasets/" + dataset_name) # this throws an error if the dataset already exists
    os.mkdir("/MAdata/sv_datasets/" + dataset_name + "/reads")
    os.mkdir("/MAdata/sv_datasets/" + dataset_name + "/alignments")
    os.mkdir("/MAdata/sv_datasets/" + dataset_name + "/calls")

    json_info_file = {
        "reference_path": reference_path,
        "create_svs_func": create_svs_func[0].__name__,
        "create_reads_funcs": [],
        "version": 1
    }

    print("creating SV-DB...")
    start = time.time()

    # create the sv_db
    database = SV_DB("/MAdata/sv_datasets/" + dataset_name + "/svs.db", "create")
    ref_pack = Pack()
    ref_pack.load(reference_path + "/ma/genome")

    # create the svs
    caller_id = create_svs_func[0](ref_pack, database, json_info_file, *create_svs_func[1])
    database.create_caller_indices()
    print(time.time() - start, "seconds")

    print("creating sequenced genome...")
    start = time.time()

    # save the sequenced genome
    seq_pack = database.reconstruct_sequenced_genome(ref_pack, caller_id)
    seq_pack.store("/MAdata/sv_datasets/" + dataset_name + "/sequenced_genome")
    with open("/MAdata/sv_datasets/" + dataset_name + "/sequenced_genome.fasta", "w") as fasta_out:
        for name, sequence in zip(seq_pack.contigNames(), seq_pack.contigSeqs()):
            fasta_out.write(">")
            fasta_out.write(name)
            fasta_out.write("\n")
            for line in textwrap.wrap(sequence, 50):
                fasta_out.write(line)
                fasta_out.write("\n")
    print(time.time() - start, "seconds")

    print("creating reads...")
    start = time.time()

    # create the reads
    for create_reads_func, name, create_reads_args in create_reads_funcs:
        for coverage in coverages:
            name_c = name + "-" + str(coverage) + "x"
            print(name_c, "...")
            json_info_file_sub = {
                "func_name": create_reads_func.__name__,
                "name": name_c,
                "coverage": coverage
            }
            create_reads_func(seq_pack, "/MAdata/sv_datasets/" + dataset_name + "/sequenced_genome.fasta", 
                            database, "/MAdata/sv_datasets/" + dataset_name + "/reads/", json_info_file_sub,
                            coverage, name_c, *create_reads_args)
            json_info_file["create_reads_funcs"].append(json_info_file_sub)
    print(time.time() - start, "seconds")

    # save the info.json file
    with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "w") as json_out:
        json.dump(json_info_file, json_out)

    print("done creating dataset:")
    print(json_info_file)


if __name__ == "__main__":
    survivor_error_profile_pac_b = "~/workspace/SURVIVOR/HG002_Pac_error_profile_bwa.txt"
    survivor_error_profile_ont = "~/workspace/SURVIVOR/NA12878_nano_error_profile_bwa.txt"

    create_dataset("/MAdata/genome/random",
                        "100nt-del-illumina-250nt-25x",
                        ( create_separate_svs, ( (sv_deletion, tuple()), 250, 1000 ) ),
                        [(create_illumina_reads_dwgsim, "illumina-250nt", (250,)),
                         (create_reads_survivor, "pacBio", (survivor_error_profile_pac_b, "pb"))],
                        [5, 10, 25, 50])
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