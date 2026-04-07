import os
import sys
import json
import re
import regex
import argparse
import subprocess
import multiprocessing as mp
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
from itertools import islice
import gzip

trans = str.maketrans("ATCGN", "TAGCN")
def rev_comp(seq): 
    return seq.translate(trans)[::-1]

def get_avg_q(qual_str):
    if not qual_str: 
        return 0
    return sum(ord(c) - 33 for c in qual_str) / len(qual_str)

def fast_dist_le2(s1, s2):
    if s1 == s2: return 0
    if abs(len(s1) - len(s2)) > 2: return 3
    return sum(1 for c1, c2 in zip(s1, s2) if c1 != c2)

def fast_lev_dist(s1, s2, max_dist=None):
    if len(s1) < len(s2):
        return fast_lev_dist(s2, s1, max_dist)
    if len(s2) == 0:
        return len(s1)
        
    if max_dist is not None and abs(len(s1) - len(s2)) > max_dist:
        return 9999
        
    prev = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        curr = [i + 1]
        min_curr = i + 1
        for j, c2 in enumerate(s2):
            ins = prev[j + 1] + 1
            dels = curr[j] + 1
            subs = prev[j] + (c1 != c2)
            val = min(ins, dels, subs)
            curr.append(val)
            if val < min_curr:
                min_curr = val
        
        if max_dist is not None and min_curr > max_dist:
            return 9999
        prev = curr
    return prev[-1]

worker_patterns = {}

def _init_pe150():
    global worker_patterns
    worker_patterns = {
        "BC1": regex.compile(r'(?:CAGGGTACGCTGTCGAGT){s<=2}(?P<bc>[ATCGN]{15,50})(?:CAGAATTCCTGCACTACG){s<=2}'),
        "BC2": regex.compile(r'(?:TGTCCAATCCATGGTGGCACT){s<=2}(?P<bc>[ATCGN]{15,50})(?:AGGATCGATATCTCGAGTCG){s<=2}'),
        "BC3": regex.compile(r'(?:GTATCCATCTTCCCACCATG){s<=2}(?P<bc>[ATCGN]{15,50})(?:CAACGACAGTACAACTACCT){s<=2}')
    }

def _init_pe300():
    global worker_patterns
    worker_patterns = {
        "BC3_PE300": regex.compile(r'(?:CTGTATCCATCTTCCCACCATG){s<=2}(?P<bc>[ATCGN]{15,50})(?:CAACGACAGTACAACTACCT){s<=2}')
    }

def _extract_bc_q(seq, qual, target):
    m = worker_patterns[target].search(seq)
    if m:
        s, e = m.span("bc")
        return m.group("bc"), get_avg_q(qual[s:e])
    rc_seq = rev_comp(seq)
    m = worker_patterns[target].search(rc_seq)
    if m:
        s, e = m.span("bc")
        rc_qual = qual[::-1]
        return m.group("bc"), get_avg_q(rc_qual[s:e])
    return None, 0

def _get_best_bc(s1, q1, s2, q2, target):
    min_q = 20
    bc1, q_v1 = _extract_bc_q(s1, q1, target)
    bc2, q_v2 = _extract_bc_q(s2, q2, target)
    if bc1 and bc2:
        if bc1 == bc2: 
            return bc1 if q_v1 >= min_q else None
        if q_v1 >= min_q and q_v1 > q_v2 + 5: 
            return bc1
        if q_v2 >= min_q and q_v2 > q_v1 + 5: 
            return bc2
        return None
    if bc1 and q_v1 >= min_q: 
        return bc1
    if bc2 and q_v2 >= min_q: 
        return bc2
    return None

def _process_pe150(args):
    chunk, t_a, t_b, p_t = args
    data_map = defaultdict(Counter)
    for s1, q1, s2, q2 in chunk:
        b_a = _get_best_bc(s1, q1, s2, q2, t_a)
        b_b = _get_best_bc(s1, q1, s2, q2, t_b)
        if b_a and b_b:
            if t_b == p_t: 
                data_map[b_b][b_a] += 1
            else: 
                data_map[b_a][b_b] += 1
    return data_map

def _process_pe300(args):
    chunk, valid_set = args
    local_map = {}
    id_map = {}
    hit_count = 0
    pattern = worker_patterns["BC3_PE300"]
    for seq, head in chunk:
        m_fwd = pattern.search(seq)
        if m_fwd:
            found_bc = m_fwd.group("bc")
            final_seq = seq
        else:
            rc_seq = rev_comp(seq)
            m_rev = pattern.search(rc_seq)
            if m_rev:
                found_bc = m_rev.group("bc")
                final_seq = rc_seq
            else:
                continue

        if found_bc in valid_set:
            if found_bc not in local_map:
                local_map[found_bc] = Counter()
                id_map[found_bc] = {}
            # 隔离 Read ID 进行纯粹的序列频次累加
            local_map[found_bc][final_seq] += 1
            if final_seq not in id_map[found_bc]:
                id_map[found_bc][final_seq] = head
            hit_count += 1
    return local_map, id_map, hit_count, len(chunk)

class scVDJSymmetricPipeline:
    def __init__(self, args):
        self.args = args
        self.sample = args.sample
        self.species = args.species
        self.out_dir = args.out_dir
        os.makedirs(self.out_dir, exist_ok=True)
        
        self.index_json = os.path.join(self.out_dir, f"{self.sample}_Clean_BC_Index.json")
        self.quant_csv = os.path.join(self.out_dir, f"{self.sample}_PE300_Quantified.csv")
        self.fs_report = os.path.join(self.out_dir, f"{self.sample}_Frameshift_Diagnosis_List.csv")
        
        self.l1_csv = os.path.join(self.out_dir, f"{self.sample}_Clone_L1_Physical.csv")
        self.l2_csv = os.path.join(self.out_dir, f"{self.sample}_Clone_L2_Clonotype.csv")
        self.l3_csv = os.path.join(self.out_dir, f"{self.sample}_Clone_L3_Lineage.csv")
        self.client_master_l4_csv = os.path.join(self.out_dir, f"{self.sample}_Client_Master_L4_Clonotypes.csv")
        self.client_top100_csv = os.path.join(self.out_dir, f"{self.sample}_Client_Top100_L4_Clonotypes.csv")
        self.report_md = os.path.join(self.out_dir, f"{self.sample}_Pipeline_Report.md")
        self.plot_dir = os.path.join(self.out_dir, "Visualizations")
        os.makedirs(self.plot_dir, exist_ok=True)
        
        self.stats = {}
        self.cores = max(4, mp.cpu_count() - 2)

    def check_files(self):
        input_files = [
            self.args.bc1_r1, self.args.bc1_r2,
            self.args.bc2_r1, self.args.bc2_r2,
            self.args.vh_fq, self.args.vk_fq
        ]
        for f in input_files:
            if not os.path.exists(f):
                print(f"[错误] 底层物理文件缺失或路径错误: {f}", flush=True)
                sys.exit(1)

    def _smart_open(self, filepath):
        with open(filepath, 'rb') as f:
            magic_num = f.read(2)
        if magic_num == b'\x1f\x8b':
            return gzip.open(filepath, 'rt')
        else:
            return open(filepath, 'rt')

    def _pe150_gen(self, r1, r2, ta, tb, pt, chunk_size=20000):
        with self._smart_open(r1) as f1, self._smart_open(r2) as f2:
            while True:
                l1 = list(islice(f1, chunk_size * 4))
                l2 = list(islice(f2, chunk_size * 4))
                if not l1 or not l2: 
                    break
                chunk = []
                for i in range(0, len(l1), 4):
                    if i+3 < len(l1) and i+3 < len(l2):
                        s1 = l1[i+1].strip().upper()
                        q1 = l1[i+3].strip()
                        s2 = l2[i+1].strip().upper()
                        q2 = l2[i+3].strip()
                        chunk.append((s1, q1, s2, q2))
                if chunk: 
                    yield (chunk, ta, tb, pt)

    def _build_pe150_map(self, r1, r2, targets, pt, name):
        print(f"[{name}] 开始多进程提取 PE150 映射...", flush=True)
        pool = mp.Pool(processes=self.cores, initializer=_init_pe150)
        gen = self._pe150_gen(r1, r2, targets[0], targets[1], pt)
        final_map = defaultdict(Counter)
        
        total_processed = 0
        for res in pool.imap_unordered(_process_pe150, gen):
            for k, v in res.items(): 
                final_map[k].update(v)
            total_processed += 1
            if total_processed % 50 == 0:
                print(f"  已处理内存块: {total_processed:,}", end="\r", flush=True)
                
        pool.close()
        pool.join()
        print(f"\n[{name}] 提取结束。共计扫描了 {total_processed:,} 个内存块。")
        return final_map

    def _resolve_bc(self, bc_map, prefix):
        f_map = {}
        stats = {f"{prefix}_Single": 0, f"{prefix}_PCR": 0, f"{prefix}_Overload": 0}
        for p_bc, tgts in bc_map.items():
            st = tgts.most_common()
            if not st: 
                continue
            top_bc, top_cnt = st[0]
            
            is_overloaded = False
            has_pcr_err = False
            
            for seq, count in st[1:]:
                if fast_dist_le2(top_bc, seq) <= 2:
                    has_pcr_err = True
                else:
                    if count / top_cnt > 0.1:
                        is_overloaded = True
            
            # 按分子(UMI)级别计数，而非原始测序深度
            if is_overloaded:
                stats[f"{prefix}_Overload"] += 1
            elif has_pcr_err:
                stats[f"{prefix}_PCR"] += 1
            else:
                stats[f"{prefix}_Single"] += 1
            
            if top_cnt >= 3 and not is_overloaded:
                if len(st) == 1 or st[1][1] / top_cnt <= 0.1:
                    f_map[p_bc] = top_bc
        return f_map, stats

    def stage1_index(self):
        print("\n[Stage 1] PE150 双轨索引构建", flush=True)
        map1 = self._build_pe150_map(self.args.bc1_r1, self.args.bc1_r2, ["BC1", "BC3"], "BC3", "BC1-BC3")
        cln1, stats1 = self._resolve_bc(map1, "BC1")
        
        map2 = self._build_pe150_map(self.args.bc2_r1, self.args.bc2_r2, ["BC2", "BC3"], "BC3", "BC3-BC2")
        cln2, stats2 = self._resolve_bc(map2, "BC2")
        
        com = set(cln1.keys()) & set(cln2.keys())
        idx_db = {b: {"BC1": cln1[b], "BC2": cln2[b]} for b in com}
        with open(self.index_json, "w") as f: 
            json.dump(idx_db, f)
            
        self.stats.update(stats1)
        self.stats.update(stats2)
        print(f"有效 BC3 枢纽总数: {len(idx_db):,}")

    def _pe300_gen(self, fq, v_set, chunk_size=50000):
        with self._smart_open(fq) as f:
            while True:
                lines = list(islice(f, chunk_size * 4))
                if not lines: 
                    break
                chunk = []
                for i in range(0, len(lines), 4):
                    if i+1 < len(lines):
                        head = lines[i].strip()
                        seq = lines[i+1].strip().upper()
                        clean_id = head[1:].split()[0].split('/')[0]
                        chunk.append((seq, clean_id))
                if chunk: 
                    yield (chunk, v_set)

    def stage2_assembly(self):
        print("\n[Stage 2] PE300 对称拓扑提取与物理组装", flush=True)
        with open(self.index_json, "r") as f: 
            idx_db = json.load(f)
        valid_set = set(idx_db.keys())

        def _scan(fq, name):
            print(f"  -> 多核并发扫描 {name}...", flush=True)
            pool = mp.Pool(processes=self.cores, initializer=_init_pe300)
            gen = self._pe300_gen(fq, valid_set)
            g_map = {k: Counter() for k in valid_set}
            g_id = {k: {} for k in valid_set}
            tot, hit = 0, 0
            
            blocks = 0
            for l_map, i_map, c_hit, c_tot in pool.imap_unordered(_process_pe300, gen):
                tot += c_tot
                hit += c_hit
                for bc, cnt in l_map.items(): 
                    g_map[bc].update(cnt)
                for bc, im in i_map.items():
                    g_id[bc].update(im)
                blocks += 1
                if blocks % 20 == 0:
                    print(f"     已扫描块: {blocks:,} (约 {tot:,} 条 Reads)", end="\r", flush=True)
                    
            pool.close()
            pool.join()
            print(f"\n     扫描完成。合计处理 {tot:,} 条序列。")
            
            best = {}
            for bc, cnt in g_map.items():
                if cnt:
                    best_seq, best_count = cnt.most_common(1)[0]
                    best_id = g_id[bc][best_seq]
                    best[bc] = {"seq": best_seq, "count": best_count, "id": best_id}
            return best, tot, hit

        vk_data, vk_t, vk_h = _scan(self.args.vk_fq, "VK_BC3")
        vh_data, vh_t, vh_h = _scan(self.args.vh_fq, "VH_BC3")

        success = 0
        with open(self.quant_csv, "w") as out:
            out.write("BC3,BC1,BC2,VH_Seq,VK_Seq,VH_Reads,VK_Reads,VH_ReadID,VK_ReadID\n")
            for bc3, info in idx_db.items():
                if bc3 in vh_data and bc3 in vk_data:
                    vh, vk = vh_data[bc3], vk_data[bc3]
                    out.write(f"{bc3},{info['BC1']},{info['BC2']},{vh['seq']},{vk['seq']},{vh['count']},{vk['count']},{vh['id']},{vk['id']}\n")
                    success += 1
        
        self.stats.update({"Raw_VH": vh_t, "Hit_VH": vh_h, "Raw_VK": vk_t, "Hit_VK": vk_h, "Total_UMIs": success})
        print(f"双链完美匹配 UMIs: {success:,} 个")

    def stage3_annotation(self):
        print("\n[Stage 3] MiXCR 解算与分型聚类", flush=True)
        df = pd.read_csv(self.quant_csv)
        
        client_features = [
            "targetSequences", "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore",
            "aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3", "aaSeqCDR3", "aaSeqFR4"
        ]

        def _run_mixcr(seqs, prefix, chain):
            fasta = os.path.join(self.out_dir, f"{prefix}_tmp.fasta")
            vdjca = os.path.join(self.out_dir, f"{prefix}_tmp.vdjca")
            anno = os.path.join(self.out_dir, f"{prefix}_tmp.tsv")
            with open(fasta, "w") as f:
                for i, s in enumerate(seqs): 
                    f.write(f">s_{i}\n{s}\n")
            
            try:
                subprocess.run(["mixcr", "-Xmx50g", "align", "-f", "-p", "generic-amplicon", "--rna", "-s", self.species, 
                                "-OsaveOriginalReads=true", "-OallowPartialAlignments=true",
                                "-OmergerParameters.minimalOverlap=10", "-OmergerParameters.minimalIdentity=0.9",
                                "--assemble-clonotypes-by", "CDR3",
                                "--floating-left-alignment-boundary", "--rigid-right-alignment-boundary", "C", 
                                fasta, vdjca], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
                               
                export_cmd = ["mixcr", "exportAlignments", "-f", "-c", chain,
                              "-readId", "-targetSequences", "-vHits", "-dHits", "-jHits",
                              "-aaFeature", "FR1", "-aaFeature", "CDR1", "-aaFeature", "FR2",
                              "-aaFeature", "CDR2", "-aaFeature", "FR3", "-aaFeature", "CDR3",
                              "-aaFeature", "FR4", "-nFeature", "CDR3", vdjca, anno]
                subprocess.run(export_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                print(f"[错误] MiXCR崩溃:\n{e.stderr.decode('utf-8', errors='ignore')}")
                sys.exit(1)
            
            anno_dict = {}
            with open(anno, "r") as f:
                header = f.readline().strip('\n').split('\t')
                def get_idx(keywords):
                    for h_idx, h_name in enumerate(header):
                        h_clean = h_name.lower().replace(" ", "").replace(".", "")
                        for kw in keywords:
                            if kw.lower() in h_clean: 
                                return h_idx
                    return -1
                
                idx_r = get_idx(["readid"])
                idx_v = get_idx(["allvhitswithscore", "vhit"])
                idx_j = get_idx(["alljhitswithscore", "jhit"])
                idx_cnt = get_idx(["nfeaturecdr3", "nseqcdr3"])
                idx_caa = get_idx(["aafeaturecdr3", "aaseqcdr3"])

                c_map = {}
                for feat in client_features:
                    if "target" in feat.lower(): 
                        c_map[feat] = get_idx(["target"])
                    elif "vhit" in feat.lower(): 
                        c_map[feat] = get_idx(["vhit", "allvhitswithscore"])
                    elif "dhit" in feat.lower(): 
                        c_map[feat] = get_idx(["dhit", "alldhitswithscore"])
                    elif "jhit" in feat.lower(): 
                        c_map[feat] = get_idx(["jhit", "alljhitswithscore"])
                    elif "aaseq" in feat.lower():
                        r_str = "FR" + feat[-1] if "fr" in feat.lower() else "CDR" + feat[-1]
                        c_map[feat] = next((i for i, c in enumerate(header) if r_str in c and ("aafeature" in c.lower() or "aaseq" in c.lower())), -1)
                    else: 
                        c_map[feat] = -1

                for line in f:
                    p = line.strip('\n').split('\t')
                    if idx_r == -1 or idx_r >= len(p): 
                        continue
                    try: 
                        rid = int(p[idx_r])
                    except ValueError: 
                        continue
                    
                    orig = seqs[rid]
                    c_nt = p[idx_cnt].upper() if idx_cnt != -1 and idx_cnt < len(p) else ""
                    c_aa = p[idx_caa] if idx_caa != -1 and idx_caa < len(p) else ""
                    
                    if c_aa and "REGION_NOT_COVERED" not in c_nt and "NOT_FOUND" not in c_nt:
                        f_dict = {feat: (p[idx] if idx != -1 and idx < len(p) else "") for feat, idx in c_map.items()}
                        anno_dict[orig] = {
                            "V_clean": p[idx_v].split('*')[0] if idx_v != -1 and idx_v < len(p) else "",
                            "J_clean": p[idx_j].split('*')[0] if idx_j != -1 and idx_j < len(p) else "",
                            "CDR3_nt": c_nt, 
                            "CDR3_aa": c_aa, 
                            **f_dict
                        }
            for file_to_del in [fasta, vdjca, anno]: 
                if os.path.exists(file_to_del): 
                    os.remove(file_to_del)
            return anno_dict

        vh_dict = _run_mixcr(df['VH_Seq'].dropna().unique().tolist(), "VH", "IGH")
        vk_dict = _run_mixcr(df['VK_Seq'].dropna().unique().tolist(), "VK", "IGK,IGL")

        for ch, d in [('VH', vh_dict), ('VK', vk_dict)]:
            df[f'{ch}_V'] = df[f'{ch}_Seq'].apply(lambda x: d.get(x, {}).get("V_clean", ""))
            df[f'{ch}_J'] = df[f'{ch}_Seq'].apply(lambda x: d.get(x, {}).get("J_clean", ""))
            df[f'{ch}_CDR3_nt'] = df[f'{ch}_Seq'].apply(lambda x: d.get(x, {}).get("CDR3_nt", ""))
            df[f'{ch}_CDR3_aa'] = df[f'{ch}_Seq'].apply(lambda x: d.get(x, {}).get("CDR3_aa", ""))
            sfx = "_H" if ch == "VH" else "_K"
            for feat in client_features:
                df[f'{feat}{sfx}'] = df[f'{ch}_Seq'].apply(lambda x: d.get(x, {}).get(feat, ""))

        df_anno = df[(df['VH_CDR3_aa'].fillna("") != "") & (df['VK_CDR3_aa'].fillna("") != "")].copy()
        if df_anno.empty:
            print("[错误] 未提取到有效双链CDR3。", flush=True)
            sys.exit(1)

        regions = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        global_fs_mask = pd.Series(False, index=df_anno.index)
        global_stop_mask = pd.Series(False, index=df_anno.index)
        
        for chain, sfx in {"VH": "_H", "VK": "_K"}.items():
            for reg in regions:
                col_name = f"aaSeq{reg}{sfx}"
                if col_name in df_anno.columns:
                    seqs = df_anno[col_name].fillna("").astype(str)
                    trunc_mask = (seqs == "") | (seqs == "region_not_covered") | (seqs.str.replace('_', '') == '')
                    clean_seqs = seqs.str.strip('_')
                    fs_mask = clean_seqs.str.contains('_', regex=False) & (~trunc_mask)
                    stop_mask = seqs.str.contains(r'\*', regex=True) & (~trunc_mask)
                    
                    self.stats[f"{chain}_{reg}_Trunc"] = trunc_mask.sum()
                    self.stats[f"{chain}_{reg}_FS"] = fs_mask.sum()
                    self.stats[f"{chain}_{reg}_Stop"] = stop_mask.sum()
                    
                    global_fs_mask = global_fs_mask | fs_mask
                    global_stop_mask = global_stop_mask | stop_mask

        self.stats["Global_FS"] = global_fs_mask.sum()
        self.stats["Global_Stop"] = global_stop_mask.sum()

        cdr3_fs = df_anno['VH_CDR3_aa'].fillna("").str.contains('_', regex=False) | df_anno['VK_CDR3_aa'].fillna("").str.contains('_', regex=False)
        cdr3_stop = df_anno['VH_CDR3_aa'].fillna("").str.contains(r'\*', regex=True) | df_anno['VK_CDR3_aa'].fillna("").str.contains(r'\*', regex=True)
        
        df_audit = df_anno.copy()
        df_audit['Defect_CDR3_Frameshift'] = cdr3_fs
        df_audit['Defect_CDR3_StopCodon'] = cdr3_stop
        df_audit['Defect_Global_Internal_FS'] = global_fs_mask
        df_audit['Defect_Global_StopCodon'] = global_stop_mask
        df_audit['Is_Productive_Functional'] = ~(cdr3_stop | cdr3_fs)
        
        audit_csv = os.path.join(self.out_dir, f"{self.sample}_Unfiltered_Library_Audit.csv")
        df_audit.to_csv(audit_csv, index=False)
        print(f"  -> 真实文库全量审计表已落盘: {os.path.basename(audit_csv)}", flush=True)

        self.stats["CDR3_FS"] = cdr3_fs.sum()
        self.stats["CDR3_Stop"] = cdr3_stop.sum()
        
        df_func = df_anno[~(cdr3_stop | cdr3_fs)].copy()
        
        if df_func.empty:
            print("[错误] 依据 CDR3 质控剔除后，无功能性序列残留。", flush=True)
            sys.exit(1)
            
        df_anno[cdr3_fs].to_csv(self.fs_report, index=False)
        
        df_func['pair_id'] = [f"Pair_{str(i+1).zfill(7)}" for i in range(len(df_func))]
        df_func['Total_Reads'] = df_func['VH_Reads'] + df_func['VK_Reads']
        
        l1_cols = ['VH_Seq', 'VK_Seq']
        df_func.groupby(l1_cols).agg(UMI_Count=('BC3', 'nunique'), Total_Reads=('Total_Reads', 'sum')).reset_index().sort_values('UMI_Count', ascending=False).to_csv(self.l1_csv, index=False)

        l2_cols = ['VH_V', 'VH_J', 'VH_CDR3_nt', 'VK_V', 'VK_J', 'VK_CDR3_nt']
        df_func.groupby(l2_cols).agg(UMI_Count=('BC3', 'nunique'), Total_Reads=('Total_Reads', 'sum')).reset_index().sort_values('UMI_Count', ascending=False).to_csv(self.l2_csv, index=False)

        l3_cols = ['VH_V', 'VH_J', 'VH_CDR3_aa', 'VK_V', 'VK_J', 'VK_CDR3_aa']
        clone_l3 = df_func.groupby(l3_cols).agg(UMI_Count=('BC3', 'nunique'), Total_Reads=('Total_Reads', 'sum')).reset_index().sort_values('UMI_Count', ascending=False)
        clone_l3.to_csv(self.l3_csv, index=False)

        print("  -> L4 广义聚类启动。基于 DP Pruning 算法进行极速降维聚合...", flush=True)
        
        df_func['VH_D'] = df_func['allDHitsWithScore_H'].fillna("").astype(str).str.split('*').str[0]
        
        l4_seed_keys = ['VH_V', 'VH_D', 'VH_J', 'VK_V', 'VK_J', 'VH_CDR3_aa', 'VK_CDR3_aa']
        idx_max = df_func.groupby(l4_seed_keys)['Total_Reads'].idxmax()
        unique_df = df_func.loc[idx_max][l4_seed_keys + ['Total_Reads']].copy()
        unique_df['rep_idx'] = idx_max.values
        
        vdj_keys = ['VH_V', 'VH_D', 'VH_J', 'VK_V', 'VK_J']
        l3_to_l4_seed = {}
        
        for _, grp in unique_df.groupby(vdj_keys):
            grp = grp.sort_values('Total_Reads', ascending=False)
            seeds = []
            for _, row in grp.iterrows():
                comb_aa = str(row['VH_CDR3_aa']) + str(row['VK_CDR3_aa'])
                matched_seed_idx = None
                
                for seed in seeds:
                    seed_aa = seed['aa']
                    max_len = max(len(comb_aa), len(seed_aa))
                    if max_len == 0:
                        continue
                        
                    allowed_dist = int(max_len * 0.1)
                    dist = fast_lev_dist(comb_aa, seed_aa, max_dist=allowed_dist)
                    
                    if dist <= allowed_dist:
                        matched_seed_idx = seed['rep_idx']
                        break
                
                if matched_seed_idx is not None:
                    l3_to_l4_seed[row['rep_idx']] = matched_seed_idx
                else:
                    seeds.append({'aa': comb_aa, 'rep_idx': row['rep_idx']})
                    l3_to_l4_seed[row['rep_idx']] = row['rep_idx']
                    
        df_func_mapped = df_func.merge(unique_df[l4_seed_keys + ['rep_idx']], on=l4_seed_keys, how='left')
        df_func_mapped['L4_Seed_Idx'] = df_func_mapped['rep_idx'].map(l3_to_l4_seed)
        
        l4_records = []
        grouped_l4 = df_func_mapped.groupby('L4_Seed_Idx')
        
        for seed_idx, grp_df in grouped_l4:
            rep_row = df_func.loc[seed_idx]
            
            bc1_set = set(grp_df['BC1'].dropna().astype(str))
            bc2_set = set(grp_df['BC2'].dropna().astype(str))
            bc3_set = set(grp_df['BC3'].dropna().astype(str))
            
            record = {
                'VH_V': rep_row['VH_V'],
                'VH_J': rep_row['VH_J'],
                'VH_CDR3_aa': rep_row['VH_CDR3_aa'],
                'VK_V': rep_row['VK_V'],
                'VK_J': rep_row['VK_J'],
                'VK_CDR3_aa': rep_row['VK_CDR3_aa'],
                'UMI_Count': len(bc3_set),
                'Total_Reads': grp_df['Total_Reads'].sum(),
                'Droplet_Count': len(bc1_set),
                'Rep_VH_Seq': rep_row['VH_Seq'],
                'Rep_VK_Seq': rep_row['VK_Seq'],
                'Rep_VH_CDR3_nt': rep_row['VH_CDR3_nt'],
                'Rep_VK_CDR3_nt': rep_row['VK_CDR3_nt'],
                'BC1': ";".join(bc1_set), 'BC2': ";".join(bc2_set), 'BC3': ";".join(bc3_set)
            }
            
            for feat in client_features:
                record[f'{feat}_H'] = rep_row.get(f'{feat}_H', '')
                record[f'{feat}_K'] = rep_row.get(f'{feat}_K', '')
                
            l4_records.append(record)
            
        df_l4 = pd.DataFrame(l4_records)
        df_l4 = df_l4.sort_values(['UMI_Count', 'Total_Reads'], ascending=[False, False]).reset_index(drop=True)
        df_l4['Clone_ID_L4'] = [f"L4_Clone_{str(i+1).zfill(6)}" for i in range(len(df_l4))]
        df_l4['Clonal_Fraction(%)'] = (df_l4['UMI_Count'] / df_l4['UMI_Count'].sum() * 100).round(4)
        df_l4['PCR_Bias'] = df_l4['Total_Reads'] / df_l4['UMI_Count']
        
        def get_defect_regions(row, chain, defect_type):
            regions = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
            defects = []
            for reg in regions:
                val = str(row.get(f'aaSeq{reg}_{chain}', ''))
                if defect_type == 'truncation':
                    if val == "" or val == "region_not_covered" or val.replace('_', '') == '':
                        defects.append(reg)
                elif defect_type == 'frameshift':
                    if '_' in val.strip('_'):
                        defects.append(reg)
                elif defect_type == 'stop':
                    if '*' in val:
                        defects.append(reg)
            return ";".join(defects) if defects else "None"

        df_l4['VH_Truncated_Regions'] = df_l4.apply(lambda r: get_defect_regions(r, 'H', 'truncation'), axis=1)
        df_l4['VK_Truncated_Regions'] = df_l4.apply(lambda r: get_defect_regions(r, 'K', 'truncation'), axis=1)
        df_l4['VH_Internal_FS_Regions'] = df_l4.apply(lambda r: get_defect_regions(r, 'H', 'frameshift'), axis=1)
        df_l4['VK_Internal_FS_Regions'] = df_l4.apply(lambda r: get_defect_regions(r, 'K', 'frameshift'), axis=1)
        df_l4['VH_StopCodon_Regions'] = df_l4.apply(lambda r: get_defect_regions(r, 'H', 'stop'), axis=1)
        df_l4['VK_StopCodon_Regions'] = df_l4.apply(lambda r: get_defect_regions(r, 'K', 'stop'), axis=1)
        
        l4_out_cols = [
            'Clone_ID_L4', 'VH_V', 'VH_J', 'VH_CDR3_aa', 'VK_V', 'VK_J', 'VK_CDR3_aa',
            'UMI_Count', 'Total_Reads', 'Droplet_Count', 
            'Rep_VH_Seq', 'Rep_VK_Seq', 'Rep_VH_CDR3_nt', 'Rep_VK_CDR3_nt',
            'Clonal_Fraction(%)', 'PCR_Bias', 
            'VH_Truncated_Regions', 'VK_Truncated_Regions',
            'VH_Internal_FS_Regions', 'VK_Internal_FS_Regions',
            'VH_StopCodon_Regions', 'VK_StopCodon_Regions',
            'targetSequences_H', 'allVHitsWithScore_H', 'allDHitsWithScore_H', 'allJHitsWithScore_H',
            'aaSeqFR1_H', 'aaSeqCDR1_H', 'aaSeqFR2_H', 'aaSeqCDR2_H', 'aaSeqFR3_H', 'aaSeqCDR3_H', 'aaSeqFR4_H',
            'targetSequences_K', 'allVHitsWithScore_K', 'allDHitsWithScore_K', 'allJHitsWithScore_K',
            'aaSeqFR1_K', 'aaSeqCDR1_K', 'aaSeqFR2_K', 'aaSeqCDR2_K', 'aaSeqFR3_K', 'aaSeqCDR3_K', 'aaSeqFR4_K'
        ]
        
        final_l4_cols = [c for c in l4_out_cols if c in df_l4.columns]
        
        df_l4[final_l4_cols].to_csv(self.client_master_l4_csv, index=False)
        print(f"  -> L4 终极总表已落盘: {os.path.basename(self.client_master_l4_csv)}", flush=True)
        
        df_l4[final_l4_cols].head(100).to_csv(self.client_top100_csv, index=False)
        print(f"  -> Top 100 克隆特征表已落盘: {os.path.basename(self.client_top100_csv)}", flush=True)
        
        self.stats.update({
            "Extracted_CDR3": len(df_anno), 
            "Final_Func": len(df_func),
            "L3_Clones": len(clone_l3), "L3_Reads": clone_l3['Total_Reads'].sum(),
            "L4_Clones": len(df_l4),
            "Top1_Bias": df_l4.iloc[0]['PCR_Bias'] if not df_l4.empty else 0,
            "Max_Bias": df_l4['PCR_Bias'].max() if not df_l4.empty else 0,
            "Mean_Bias": df_l4['Total_Reads'].sum() / df_l4['UMI_Count'].sum() if not df_l4.empty else 0
        })
        print(f"[Stage 3] 发现 {len(clone_l3):,} 个原始 L3 克隆，聚类收敛为 {len(df_l4):,} 个 L4 终极克隆型。", flush=True)

    def plot_and_report(self):
        print("\n[Stage 4] 全景图表与质控报告生成", flush=True)
        
        if not os.path.exists(self.client_master_l4_csv):
            return

        df = pd.read_csv(self.client_master_l4_csv)
        if df.empty or df['UMI_Count'].sum() == 0:
            return
        
        counts = df['UMI_Count']
        p = counts / counts.sum()
        shannon = -np.sum(p * np.log2(p))
        sorted_counts = np.sort(counts)
        cum_counts = np.cumsum(sorted_counts)
        n = len(counts)
        
        gini = (n + 1 - 2 * np.sum(cum_counts) / cum_counts[-1]) / n
        top10_freq = p.sort_values(ascending=False).head(10).sum() * 100

        sns.set_theme(style="whitegrid", font_scale=1.2)
        fig, ax = plt.subplots(figsize=(12, 6))
        df_top = df.head(20).copy()
        df_top['Clone_ID'] = df_top['Clone_ID_L4']
        sns.barplot(data=df_top, x='Clone_ID', y='Clonal_Fraction(%)', color='steelblue', ax=ax)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
        ax.set_title("Top 20 Clonal Fractions (L4 Clusters)")
        plt.tight_layout()
        fig.savefig(os.path.join(self.plot_dir, "01_Top20_L4_Clones.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

        df_sorted = df.sort_values(by='UMI_Count', ascending=False).reset_index(drop=True)
        df_sorted['Cum'] = df_sorted['UMI_Count'].cumsum() / df_sorted['UMI_Count'].sum() * 100
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.plot(np.arange(1, len(df_sorted)+1), df_sorted['Cum'], color='#d95f02', linewidth=2.5)
        ax.set_xscale('log')
        ax.set_title("Clonal Accumulation Curve")
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, "04_Accumulation_Curve.png"), dpi=300, bbox_inches='tight')
        plt.close(fig)

        s = self.stats
        tot_bc1 = s.get('BC1_Single', 0) + s.get('BC1_PCR', 0) + s.get('BC1_Overload', 0)
        tot_bc2 = s.get('BC2_Single', 0) + s.get('BC2_PCR', 0) + s.get('BC2_Overload', 0)
        total_seqs = s.get('Extracted_CDR3', 1) or 1
        
        regions_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
        
        vh_fs_rates = [(s.get(f'VH_{r}_FS', 0) / total_seqs) * 100 for r in regions_order]
        vk_fs_rates = [(s.get(f'VK_{r}_FS', 0) / total_seqs) * 100 for r in regions_order]
        vh_stop_rates = [(s.get(f'VH_{r}_Stop', 0) / total_seqs) * 100 for r in regions_order]
        vk_stop_rates = [(s.get(f'VK_{r}_Stop', 0) / total_seqs) * 100 for r in regions_order]
        vh_trunc_rates = [(s.get(f'VH_{r}_Trunc', 0) / total_seqs) * 100 for r in regions_order]
        vk_trunc_rates = [(s.get(f'VK_{r}_Trunc', 0) / total_seqs) * 100 for r in regions_order]
        
        heatmap_matrix = np.array([
            vh_trunc_rates, vk_trunc_rates, 
            vh_fs_rates, vk_fs_rates, 
            vh_stop_rates, vk_stop_rates
        ])
        y_labels = [
            'VH Truncation', 'VK Truncation', 
            'VH Internal FS', 'VK Internal FS', 
            'VH Stop Codon', 'VK Stop Codon'
        ]
        
        sns.set_theme(style="whitegrid", font_scale=1.1)
        fig_mut, (ax_fs, ax_trunc, ax_heat) = plt.subplots(3, 1, figsize=(12, 14), sharex=True, gridspec_kw={'height_ratios': [1, 1, 1.8]})
        
        ax_fs.plot(regions_order, vh_fs_rates, marker='o', markersize=7, label='VH Internal FS', color='#d95f02', linewidth=2.5)
        ax_fs.plot(regions_order, vk_fs_rates, marker='s', markersize=7, label='VK Internal FS', color='#1b9e77', linewidth=2.5)
        ax_fs.plot(regions_order, vh_stop_rates, marker='^', markersize=7, label='VH Stop Codon', color='#d95f02', linestyle='--', alpha=0.7)
        ax_fs.plot(regions_order, vk_stop_rates, marker='v', markersize=7, label='VK Stop Codon', color='#1b9e77', linestyle='--', alpha=0.7)
        ax_fs.set_ylabel('Mutation Rate (%)')
        ax_fs.set_title('Pathology Hotspots (Internal Frameshift & Stop Codon)', pad=10)
        ax_fs.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=4, fontsize=9)
        
        ax_trunc.plot(regions_order, vh_trunc_rates, marker='o', markersize=7, label='VH Truncation (Missing)', color='#7570b3', linewidth=2.5)
        ax_trunc.plot(regions_order, vk_trunc_rates, marker='s', markersize=7, label='VK Truncation (Missing)', color='#e7298a', linewidth=2.5)
        ax_trunc.set_ylabel('Truncation Rate (%)')
        ax_trunc.set_title('Physical Truncation & Coverage Drop-off', pad=10)
        ax_trunc.legend(loc='upper right', fontsize=10)
        ax_trunc.set_ylim(-5, 105)
        
        sns.heatmap(heatmap_matrix, annot=True, fmt=".2f", cmap="YlOrRd", vmax=20,
                    xticklabels=regions_order, yticklabels=y_labels, 
                    cbar_kws={'label': 'Defect Percentage (%)'}, ax=ax_heat, linewidths=.5)
        ax_heat.set_xlabel('Antibody Topologic Regions (5\' -> 3\')', labelpad=15)
        ax_heat.set_title('Comprehensive Defect Density Matrix', pad=10)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.plot_dir, "05_Mutation_Hotspots_Profiler.png"), dpi=300, bbox_inches='tight')
        plt.close(fig_mut)
        
        regional_md = "\n### 高分辨率区域病理分布 (High-Res Regional Diagnostics)\n"
        regional_md += "> **说明**: “缺失/截断”代表由于测序边界或引物覆盖不足导致的自然截断；“内部移码”已剥离边缘截断伪影，反映真实的聚合酶滑移。\n\n"
        regional_md += "| Region | VH 缺失/截断 | VH 内部移码 (_) | VH 终止 (*) | VK 缺失/截断 | VK 内部移码 (_) | VK 终止 (*) |\n"
        regional_md += "|---|---|---|---|---|---|---|\n"
        for reg in regions_order:
            vh_trunc, vh_fs, vh_stop = s.get(f'VH_{reg}_Trunc', 0), s.get(f'VH_{reg}_FS', 0), s.get(f'VH_{reg}_Stop', 0)
            vk_trunc, vk_fs, vk_stop = s.get(f'VK_{reg}_Trunc', 0), s.get(f'VK_{reg}_FS', 0), s.get(f'VK_{reg}_Stop', 0)
            
            vh_t_p, vh_f_p, vh_s_p = (vh_trunc/total_seqs)*100, (vh_fs/total_seqs)*100, (vh_stop/total_seqs)*100
            vk_t_p, vk_f_p, vk_s_p = (vk_trunc/total_seqs)*100, (vk_fs/total_seqs)*100, (vk_stop/total_seqs)*100
            
            regional_md += f"| {reg} | {vh_trunc:,} ({vh_t_p:.2f}%) | {vh_fs:,} ({vh_f_p:.2f}%) | {vh_stop:,} ({vh_s_p:.2f}%) | {vk_trunc:,} ({vk_t_p:.2f}%) | {vk_fs:,} ({vk_f_p:.2f}%) | {vk_stop:,} ({vk_s_p:.2f}%) |\n"
        
        md = f"""
# {self.sample} 全长单细胞 VDJ 分析总决报告 (L4-Optimized)

| Pipeline Stage | 关键数据节点 | Counts | 转化率/占比 (%) | 备注 |
|---|---|---|---|---|
| Merged Input & Capture | VH 总测序深度 | {s.get('Raw_VH', 0):,} | 100.00% | 物理组装前初始数据量 |
| | VK 总测序深度 | {s.get('Raw_VK', 0):,} | 100.00% | |
| | VH 命中白名单 Reads | {s.get('Hit_VH', 0):,} | {s.get('Hit_VH', 0)/(s.get('Raw_VH', 1) or 1)*100:.2f}% | |
| | VK 命中白名单 Reads | {s.get('Hit_VK', 0):,} | {s.get('Hit_VK', 0)/(s.get('Raw_VK', 1) or 1)*100:.2f}% | |
| Collision Resolution | BC1 端纯净分子枢纽 (UMIs) | {s.get('BC1_Single', 0):,} | {s.get('BC1_Single', 0)/(tot_bc1 or 1)*100:.2f}% | |
| | BC1 端含PCR错配枢纽 (UMIs) | {s.get('BC1_PCR', 0):,} | {s.get('BC1_PCR', 0)/(tot_bc1 or 1)*100:.2f}% | |
| | BC1 端物理过载枢纽 (UMIs) | {s.get('BC1_Overload', 0):,} | {s.get('BC1_Overload', 0)/(tot_bc1 or 1)*100:.2f}% | |
| | BC2 端纯净分子枢纽 (UMIs) | {s.get('BC2_Single', 0):,} | {s.get('BC2_Single', 0)/(tot_bc2 or 1)*100:.2f}% | |
| | BC2 端含PCR错配枢纽 (UMIs) | {s.get('BC2_PCR', 0):,} | {s.get('BC2_PCR', 0)/(tot_bc2 or 1)*100:.2f}% | |
| | BC2 端物理过载枢纽 (UMIs) | {s.get('BC2_Overload', 0):,} | {s.get('BC2_Overload', 0)/(tot_bc2 or 1)*100:.2f}% | |
| Physical Assembly | UMIs | {s.get('Total_UMIs', 0):,} | 100.00% | 物理桥接成功的分子总数 |
| MiXCR Annotation | 成功提取双链 CDR3 的分子 | {s.get('Extracted_CDR3', 0):,} | {s.get('Extracted_CDR3', 0)/(s.get('Total_UMIs', 1) or 1)*100:.2f}% | |
| Global Mutation Profile | 全局存在真移码 (_) 的分子 | {s.get('Global_FS', 0):,} | {s.get('Global_FS', 0)/(s.get('Extracted_CDR3', 1) or 1)*100:.2f}% | 仅统计非边缘截断导致的内部移码 |
| | 全局存在终止 (*) 的分子 | {s.get('Global_Stop', 0):,} | {s.get('Global_Stop', 0)/(s.get('Extracted_CDR3', 1) or 1)*100:.2f}% | |
| Functional Filtering | CDR3 移码 (执行剔除) | {s.get('CDR3_FS', 0):,} | {s.get('CDR3_FS', 0)/(s.get('Extracted_CDR3', 1) or 1)*100:.2f}% | |
| | CDR3 终止 (执行剔除) | {s.get('CDR3_Stop', 0):,} | {s.get('CDR3_Stop', 0)/(s.get('Extracted_CDR3', 1) or 1)*100:.2f}% | |
| | 最终全长功能分子 | {s.get('Final_Func', 0):,} | {s.get('Final_Func', 0)/(s.get('Extracted_CDR3', 1) or 1)*100:.2f}% | |
| Biological Clustering | Level 3 进化论谱系金标准 | {s.get('L3_Clones', 0):,} | - | 基于 100% CDR3 氨基酸匹配 |
| | Level 4 终极聚合克隆型 | {s.get('L4_Clones', 0):,} | - | 基于同 V/D/J 且 >= 90% AA 匹配 |
| Diversity Profiling | Total Reads | {s.get('L3_Reads', 0):,} | - | |
| | Shannon Index | {shannon:.4f} | - | |
| | Gini Index | {gini:.4f} | - | |
| | L4 Top 10 clone freq | {top10_freq:.2f}% | - | |
| | Top 1 Clone PCR bias | {s.get('Top1_Bias', 0):.2f} | - | (Reads/UMI) |
| | Max PCR bias | {s.get('Max_Bias', 0):.2f} | - | (Reads/UMI) |
| | Mean PCR bias | {s.get('Mean_Bias', 0):.2f} | - | (Reads/UMI) |
{regional_md}
"""
        with open(self.report_md, "w", encoding="utf-8") as f: 
            f.write(md.strip())

    def run(self):
        self.check_files()
        self.stage1_index()
        self.stage2_assembly()
        self.stage3_annotation()
        self.plot_and_report()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", required=True)
    parser.add_argument("-sp", "--species", default="hsa")
    parser.add_argument("--bc1_r1", required=True)
    parser.add_argument("--bc1_r2", required=True)
    parser.add_argument("--bc2_r1", required=True)
    parser.add_argument("--bc2_r2", required=True)
    parser.add_argument("--vh_fq", required=True)
    parser.add_argument("--vk_fq", required=True)
    parser.add_argument("-o", "--out_dir", required=True)
    args = parser.parse_args()
    
    scVDJSymmetricPipeline(args).run()