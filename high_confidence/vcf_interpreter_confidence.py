from caller_analysis.vcf_interpreters import *

min_sv_size = 5

def confidence_call_interpreter(call, call_inserter, pack, error_file, call_desc):
    def find_confidence(call):
        return int(float(call["QUAL"])*100)

    try:
        from_pos = int(call["POS"])
        chrom = call["CHROM"]
        name_trans_dict = {
                "chr1": "CM000663.2",
                "chr2": "CM000664.2",
                "chr3": "CM000665.2",
                "chr4": "CM000666.2",
                "chr5": "CM000667.2",
                "chr6": "CM000668.2",
                "chr7": "CM000669.2",
                "chr8": "CM000670.2",
                "chr9": "CM000671.2",
                "chr10": "CM000672.2",
                "chr11": "CM000673.2",
                "chr12": "CM000674.2",
                "chr13": "CM000675.2",
                "chr14": "CM000676.2",
                "chr15": "CM000677.2",
                "chr16": "CM000678.2",
                "chr17": "CM000679.2",
                "chr18": "CM000680.2",
                "chr19": "CM000681.2",
                "chr20": "CM000682.2",
                "chr21": "CM000683.2",
                "chr22": "CM000684.2",
                "chrX": "CM000685.2",
                "chrY": "CM000686.2"
            }
        if chrom in name_trans_dict:
            chrom = name_trans_dict[chrom]
        from_pos += pack.start_of_sequence(chrom)
        size_ref = len(call["REF"])
        size_alt = len(call["ALT"])
        if size_ref >= min_sv_size or size_alt >= min_sv_size:
            call_to_insert = SvCall(from_pos, from_pos + size_ref, 1, 1, False, find_confidence(call), 1)
            call_inserter.insert(call_to_insert)
            call_desc.insert(call_to_insert.id, call["INFO"]["datasetnames"])

    except Exception as e:
        log_error(call, error_file, "high_confidence", e)