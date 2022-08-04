import pgrtk
import os
from collections import Counter


ref_db =pgrtk.AGCFile("pgr-tk-HGRP-y1-evaluation-set-v0.agc")
sdb = pgrtk.SeqIndexDB()
sdb.load_from_agc_index("pgr-tk-HGRP-y1-evaluation-set-v0")

seq_info = sdb.seq_info.copy()

CMRG_coordinates = {}

with open("input.bed") as f:
    for r in f:
        r = r.strip().split("\t")
        n = "chr{}_hg38".format(r[0])
        CMRG_coordinates[r[3]]=(n, int(r[1]), int(r[2]))

CMRG_hg38_seq = {}
padding = 10000
for g, c in CMRG_coordinates.items():
    print(c[0])
    seq = ref_db.get_sub_seq('hg38_tagged.fa', c[0], c[1]-padding, c[2]+padding)
    CMRG_hg38_seq[g] = (c, c[2] - c[1], seq)

out_dir = "output_dir"
os.makedirs(out_dir, exist_ok = True)
f1 = open(os.path.join(out_dir, "out_seqs_with_padding.fa"), "w") 
f2 = open(os.path.join(out_dir, "out_seqs.fa"), "w") 

for g in CMRG_hg38_seq:
    _, q_len, q_seq = CMRG_hg38_seq[g]
    aln_range = pgrtk.query_sdb(sdb, q_seq, merge_range_tol=len(q_seq) * 0.25)
    n_copy = {}
    seq_list = []
    for k in list(aln_range.keys()):
        b, e = aln_range[k][0][0:2]
        if e-b < len(q_seq) * 0.25:
            continue
        n_copy[k] = len(aln_range[k])
        copy_count = Counter(n_copy.values())
        ctg, src, len_ = sdb.seq_info[k]
        rgns = aln_range[k].copy()
        rgns = pgrtk.merge_regions(rgns, tol=int(len(q_seq)*0.25))
        for rgn in rgns:
            b, e, length, orientation, aln = rgn
            if length < len(q_seq)*0.25:
                continue
            if length > len(q_seq) * 10:
                continue
            seq =  sdb.get_sub_seq(src, ctg, b, e)
            s_name = "{}_{}_{}_{}_{} {}".format(g, ctg, b, e, orientation, src)
            print(f">{s_name}", file=f1)
            print(pgrtk.u8_to_string(seq), file=f1)
            seq_list.append((s_name, seq))

    for nc, nh in copy_count.items():
        print("#NC", g, nc, nh, sep="\t")

    new_sdb = pgrtk.SeqIndexDB() 
    new_sdb.load_from_seq_list(seq_list, w=48, k=48, r=1, min_span=24)
    new_sdb.generate_mapg_gfa(0, os.path.join(out_dir, f"{g}_48_48_1_24.gfa"))
    new_sdb.write_midx_to_text_file(os.path.join(out_dir, f"{g}_48_48_1_24.gidx"))


    interval = (padding, padding + q_len,)
    pos_map = sdb.map_positions_in_seq(interval, q_seq, 0.01, 32, 32, 32, 1000)
    #for p in pos_map:
    #    print("X", g, p, sep="\t")
    m_itvl = pgrtk.map_intervals_in_sdb(sdb, interval, q_seq, 0.01)
    print("#NI", g, interval[0], interval[1], len(m_itvl), sep="\t")

    for sid in m_itvl:
        strand, left_p, right_p = m_itvl[sid]
        ctg, src, t_len = seq_info[sid]
        print(g, src, ctg, strand, left_p, right_p, sep="\t")
        if strand == 1:
            b, e = t_len - right_p, t_len - left_p
        else:
            b, e = left_p, right_p 

        if b < e and  e - b < len(q_seq) * 10:
            if e >= t_len:
                e = t_len - 1
            if b >= t_len:
                b = t_len - 1
            if not b < e:
                continue
            seq =  sdb.get_sub_seq(src, ctg, b, e)
            if strand == 1:
                seq = pgrtk.rc_byte_seq(seq)
            print(">{}_{}_{}_{}_{} {}".format(g, ctg, b, e, strand, src), file=f2)
            print(pgrtk.u8_to_string(seq), file=f2)
f1.close()
f2.close()
