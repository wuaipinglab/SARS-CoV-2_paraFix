import os
import glob
import re
import shutil


OUT_DIR = "Data/sampled_trees_with_MSA/"

for seq_fp in glob.glob(os.path.join("Data/sampled_trees/", "*", "*.fasta")):
    dn = os.path.dirname(seq_fp)
    seq_fn = os.path.basename(seq_fp)
    (dt,) = re.findall(r"[0-9]{8}", seq_fn)
    dt = "-".join([dt[0:4], dt[4:6], dt[6:8]])
    print(seq_fn, dt)

    out_dir = os.path.join(OUT_DIR, dt)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_seq_fp = os.path.join(out_dir, dt + ".fasta")
    if not os.path.exists(out_seq_fp):
        shutil.copyfile(seq_fp, out_seq_fp)

    tree_fn = seq_fn + ".treefile"
    tree_fp = os.path.join(dn, tree_fn)
    out_tree_fp = os.path.join(out_dir, dt + ".nwk")
    if not os.path.exists(out_tree_fp):
        shutil.copyfile(tree_fp, out_tree_fp)
        
    out_seq_aa_fp = os.path.join(out_dir, dt + "_aa.fasta")
    if not os.path.exists(out_seq_aa_fp):
        shutil.copyfile(
            os.path.join(out_dir, os.path.splitext(seq_fn)[0] + "_aa.fasta"),
            out_seq_aa_fp
        )
        
# OUT_DIR = "Data/latest_trees_with_MSA/"
# 
# 
# for tree_fp in glob.glob(os.path.join(OUT_DIR, "*", "*.treefile")):
#     dn = os.path.dirname(tree_fp)
#     out_tree_fp = os.path.join(dn, os.path.basename(dn) + ".nwk")
#     shutil.copyfile(
#         tree_fp,
#         out_tree_fp
#     )
#     print(tree_fp)
