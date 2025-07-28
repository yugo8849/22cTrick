#!/usr/bin/env python3
"""
22c-trick Saturation PCR Designer - Streamlit GUI版
NEBuilder最適化プライマー設計ツール
"""

import streamlit as st
import pandas as pd
import itertools
import math
from typing import List, Tuple, Dict, Optional
import io
import base64

class OptimizedSaturationDesigner:
    def __init__(self):
        # 遺伝コード
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        # 22c-trick コドンセット
        self.codon_sets = {
            'NDT': {'codons': ['ATT', 'AAT', 'ACT', 'GGT', 'GAT', 'GTT', 'CGT', 'CAT', 'CTT', 'TTT', 'TAT', 'TGT'], 'count': 12},
            'VHG': {'codons': ['ATG', 'ACG', 'AAG', 'GAG', 'GCG', 'GTG', 'CAG', 'CCG', 'CTG'], 'count': 9},
            'TGG': {'codons': ['TGG'], 'count': 1}
        }
        
        # 相補塩基
        self.complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'D': 'H', 'V': 'B', 'H': 'D', 'B': 'V'}
        
        # デジェネレート塩基
        self.degenerate_bases = {
            'N': ['A', 'T', 'G', 'C'],
            'D': ['A', 'T', 'G'],
            'V': ['A', 'C', 'G'],
            'H': ['A', 'T', 'C'],
            'B': ['T', 'G', 'C']
        }
        
        # 一般的な制限酵素データベース
        self.restriction_enzymes = {
            'EcoRI': {'sequence': 'GAATTC', 'cut_pos': 1},
            'BamHI': {'sequence': 'GGATCC', 'cut_pos': 1}, 
            'HindIII': {'sequence': 'AAGCTT', 'cut_pos': 1},
            'XhoI': {'sequence': 'CTCGAG', 'cut_pos': 1},
            'SalI': {'sequence': 'GTCGAC', 'cut_pos': 1},
            'XbaI': {'sequence': 'TCTAGA', 'cut_pos': 1},
            'SpeI': {'sequence': 'ACTAGT', 'cut_pos': 1},
            'NotI': {'sequence': 'GCGGCCGC', 'cut_pos': 2},
            'PstI': {'sequence': 'CTGCAG', 'cut_pos': 5},
            'SacI': {'sequence': 'GAGCTC', 'cut_pos': 5},
            'KpnI': {'sequence': 'GGTACC', 'cut_pos': 5},
            'SmaI': {'sequence': 'CCCGGG', 'cut_pos': 3},
            'EcoRV': {'sequence': 'GATATC', 'cut_pos': 3},
            'HaeIII': {'sequence': 'GGCC', 'cut_pos': 2},
            'ApaI': {'sequence': 'GGGCCC', 'cut_pos': 5},
            'BglII': {'sequence': 'AGATCT', 'cut_pos': 1},
            'ClaI': {'sequence': 'ATCGAT', 'cut_pos': 2},
            'DraI': {'sequence': 'TTTAAA', 'cut_pos': 3},
            'NcoI': {'sequence': 'CCATGG', 'cut_pos': 1},
            'NdeI': {'sequence': 'CATATG', 'cut_pos': 2}
        }
    
    def find_restriction_sites(self, template_seq: str) -> Dict[str, List[int]]:
        """テンプレート配列中の制限酵素サイトを検索"""
        sites = {}
        
        for enzyme_name, enzyme_data in self.restriction_enzymes.items():
            recognition_seq = enzyme_data['sequence']
            positions = []
            
            # 順方向検索
            start = 0
            while True:
                pos = template_seq.find(recognition_seq, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1
            
            # 逆相補も検索
            rev_comp_seq = self.reverse_complement(recognition_seq)
            if rev_comp_seq != recognition_seq:  # 回文配列でない場合
                start = 0
                while True:
                    pos = template_seq.find(rev_comp_seq, start)
                    if pos == -1:
                        break
                    positions.append(pos)
                    start = pos + 1
            
            if positions:
                sites[enzyme_name] = sorted(set(positions))
        
        return sites
    
    def design_pcr_oligos(self, primers: List[Dict], gene_name: str) -> List[Dict]:
        """PCR用オリゴ設計（既存プライマーの共通配列を利用）"""
        pcr_oligos = []
        
        if not primers:
            return pcr_oligos
        
        # フォワード・リバースそれぞれで共通配列を抽出
        forward_seqs = [p['forward_sequence'] for p in primers]
        reverse_seqs = [p['reverse_sequence'] for p in primers]
        
        # 5'末端から共通配列を検索
        def find_common_prefix(sequences):
            if not sequences:
                return ""
            min_len = min(len(seq) for seq in sequences)
            for i in range(min_len):
                if not all(seq[i] == sequences[0][i] for seq in sequences):
                    return sequences[0][:i]
            return sequences[0][:min_len]
        
        # 3'末端から共通配列を検索  
        def find_common_suffix(sequences):
            if not sequences:
                return ""
            min_len = min(len(seq) for seq in sequences)
            for i in range(1, min_len + 1):
                if not all(seq[-i] == sequences[0][-i] for seq in sequences):
                    return sequences[0][-(i-1):] if i > 1 else ""
            return sequences[0][-min_len:]
        
        # フォワードプライマー用PCRオリゴ
        forward_prefix = find_common_prefix(forward_seqs)
        forward_suffix = find_common_suffix(forward_seqs)
        
        if len(forward_prefix) >= 18:  # 最小アニーリング長
            pcr_oligos.append({
                'name': f"{gene_name}_PCR_Vector_Fwd",
                'sequence': forward_prefix,
                'length': len(forward_prefix),
                'tm': round(self.calculate_tm(forward_prefix), 1),
                'type': 'PCR_Forward',
                'purpose': 'Vector amplification (Forward)',
                'description': f'Common 5\' sequence from saturation primers'
            })
        
        # フルレングス共通配列も追加（より確実なPCR用）
        if len(forward_seqs[0]) >= 25:
            pcr_oligos.append({
                'name': f"{gene_name}_PCR_Vector_Fwd_Full",
                'sequence': forward_seqs[0],  # 代表配列を使用
                'length': len(forward_seqs[0]),
                'tm': round(self.calculate_tm(forward_seqs[0]), 1),
                'type': 'PCR_Forward_Full',
                'purpose': 'Vector amplification (Forward, full)',
                'description': f'Representative forward primer sequence'
            })
        
        # リバースプライマー用PCRオリゴ
        reverse_prefix = find_common_prefix(reverse_seqs)
        
        if len(reverse_prefix) >= 18:
            pcr_oligos.append({
                'name': f"{gene_name}_PCR_Vector_Rev",
                'sequence': reverse_prefix,
                'length': len(reverse_prefix),
                'tm': round(self.calculate_tm(reverse_prefix), 1),
                'type': 'PCR_Reverse',
                'purpose': 'Vector amplification (Reverse)',
                'description': f'Common 5\' sequence from saturation primers'
            })
        
        if len(reverse_seqs[0]) >= 25:
            pcr_oligos.append({
                'name': f"{gene_name}_PCR_Vector_Rev_Full",
                'sequence': reverse_seqs[0],
                'length': len(reverse_seqs[0]),
                'tm': round(self.calculate_tm(reverse_seqs[0]), 1),
                'type': 'PCR_Reverse_Full',
                'purpose': 'Vector amplification (Reverse, full)',
                'description': f'Representative reverse primer sequence'
            })
        
        return pcr_oligos
    
    def design_restriction_oligos(self, template_seq: str, restriction_sites: Dict[str, List[int]], 
                                 selected_enzyme: str, selected_position: int, 
                                 mutation_positions: List[int], gene_name: str, 
                                 overlap_length: int = 20) -> List[Dict]:
        """制限酵素切断用のオリゴ設計"""
        
        if selected_enzyme not in self.restriction_enzymes:
            return []
        
        enzyme_data = self.restriction_enzymes[selected_enzyme]
        recognition_seq = enzyme_data['sequence']
        
        # 制限酵素サイト + 20bp上流を相同配列とする
        upstream_start = max(0, selected_position - 20)
        restriction_end = selected_position + len(recognition_seq)
        
        # 相同配列
        homology_seq = template_seq[upstream_start:restriction_end]
        
        # NEBuilder用オリゴ設計（制限酵素サイトを利用）
        restriction_oligos = []
        
        # 制限酵素サイト情報
        restriction_oligos.append({
            'name': f"{gene_name}_Restriction_Info",
            'sequence': f"{selected_enzyme}: {recognition_seq} at position {selected_position}",
            'length': len(recognition_seq),
            'tm': 'N/A',
            'type': 'Restriction_Info',
            'purpose': 'Restriction enzyme information',
            'description': f'Cut with {selected_enzyme} at position {selected_position}'
        })
        
        # ベクター用フォワードオリゴ（制限酵素サイト含む）
        restriction_oligos.append({
            'name': f"{gene_name}_Vector_Fwd_{selected_enzyme}",
            'sequence': homology_seq,
            'length': len(homology_seq),
            'tm': round(self.calculate_tm(homology_seq), 1),
            'type': 'Vector_Forward',
            'purpose': 'Vector linearization (Forward)',
            'description': f'20bp upstream + {selected_enzyme} site'
        })
        
        # 対応するリバースオリゴも設計
        reverse_homology = self.reverse_complement(homology_seq)
        restriction_oligos.append({
            'name': f"{gene_name}_Vector_Rev_{selected_enzyme}",
            'sequence': reverse_homology,
            'length': len(reverse_homology),
            'tm': round(self.calculate_tm(reverse_homology), 1),
            'type': 'Vector_Reverse',
            'purpose': 'Vector linearization (Reverse)',
            'description': f'Reverse complement of forward vector oligo'
        })
        
        return restriction_oligos
    
    def translate_dna(self, dna_seq: str) -> str:
        """DNA配列をアミノ酸配列に翻訳"""
        protein = ""
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3]
            if len(codon) == 3:
                protein += self.genetic_code.get(codon, 'X')
        return protein
        """DNA配列をアミノ酸配列に翻訳"""
        protein = ""
        for i in range(0, len(dna_seq), 3):
            codon = dna_seq[i:i+3]
            if len(codon) == 3:
                protein += self.genetic_code.get(codon, 'X')
        return protein
    
    def reverse_complement(self, seq: str) -> str:
        """逆相補配列を取得"""
        return ''.join(self.complement.get(base, base) for base in reversed(seq))
    
    def calculate_tm(self, seq: str) -> float:
        """Tm計算（Wallace rule + GC補正）"""
        # デジェネレート塩基の平均GC含量を計算
        total_gc = 0
        total_at = 0
        
        for base in seq:
            if base in self.degenerate_bases:
                possible_bases = self.degenerate_bases[base]
                gc_count = sum(1 for b in possible_bases if b in ['G', 'C'])
                at_count = sum(1 for b in possible_bases if b in ['A', 'T'])
                total_gc += gc_count / len(possible_bases)
                total_at += at_count / len(possible_bases)
            elif base in ['G', 'C']:
                total_gc += 1
            elif base in ['A', 'T']:
                total_at += 1
        
        if len(seq) < 14:
            return 2 * total_at + 4 * total_gc
        else:
            gc_content = total_gc / len(seq) if len(seq) > 0 else 0
            return 64.9 + 41 * (gc_content - 16.4 / len(seq))
    
    def design_nebuilder_primers(self, template_seq: str, mutation_positions: List[int], 
                                gene_name: str = "Gene", overlap_length: int = 20,
                                vector_method: str = "PCR", restriction_enzyme: str = None,
                                restriction_position: int = None) -> Tuple[List[Dict], List[Dict]]:
        """NEBuilder最適化プライマー設計（ベクター準備方法対応）"""
        
        # Saturation mutagenesis用プライマー設計
        if len(mutation_positions) == 1:
            sat_primers = self._design_single_nebuilder(template_seq, mutation_positions[0], gene_name, overlap_length)
        elif len(mutation_positions) == 2:
            sat_primers = self._design_double_nebuilder(template_seq, mutation_positions, gene_name, overlap_length)
        else:
            raise ValueError("1〜2箇所の変異位置のみ対応")
        
        # ベクター用オリゴ設計
        vector_oligos = []
        
        if vector_method == "PCR":
            # PCR用オリゴ設計
            vector_oligos = self.design_pcr_oligos(sat_primers, gene_name)
        elif vector_method == "Restriction" and restriction_enzyme and restriction_position is not None:
            # 制限酵素用オリゴ設計
            restriction_sites = self.find_restriction_sites(template_seq)
            vector_oligos = self.design_restriction_oligos(
                template_seq, restriction_sites, restriction_enzyme, 
                restriction_position, mutation_positions, gene_name, overlap_length
            )
        
        return sat_primers, vector_oligos
    
    def _find_optimal_annealing_region(self, template_seq: str, center_pos: int, min_length: int = 18, 
                                     min_tm: float = 55.0, max_total: int = 60) -> Tuple[int, int]:
        """最適なアニーリング領域を見つける"""
        best_start, best_end = center_pos - min_length//2, center_pos + min_length//2
        
        # 最小条件から開始して、Tm条件を満たすまで延長
        for length in range(min_length, max_total - 20):  # オーバーラップ分を引く
            start = max(0, center_pos - length//2)
            end = min(len(template_seq), start + length)
            
            if end - start < min_length:
                continue
                
            test_seq = template_seq[start:end]
            tm = self.calculate_tm(test_seq)
            
            if tm >= min_tm:
                return start, end
                
        return best_start, best_end
    
    def _design_single_nebuilder(self, template_seq: str, position: int, gene_name: str, overlap_length: int) -> List[Dict]:
        """1箇所変異用NEBuilderプライマー"""
        primers = []
        
        # アニーリング領域の最適化
        center = position + 1  # 変異位置の中央
        anneal_start, anneal_end = self._find_optimal_annealing_region(template_seq, center)
        
        primer_designs = [('NDT', 12, 1), ('VHG', 9, 2), ('TGG', 1, 3)]
        
        for codon_type, combinations, primer_num in primer_designs:
            # 変異を含むアニーリング領域を作成
            left_part = template_seq[anneal_start:position]
            right_part = template_seq[position + 3:anneal_end]
            annealing_region = left_part + codon_type + right_part
            
            # オーバーラップ領域を追加
            left_overlap = template_seq[max(0, anneal_start - overlap_length):anneal_start]
            right_overlap = template_seq[anneal_end:min(len(template_seq), anneal_end + overlap_length)]
            
            forward_seq = left_overlap + annealing_region + right_overlap
            reverse_seq = self.reverse_complement(forward_seq)
            
            # 長さ制限チェック
            if len(forward_seq) > 60:
                # 長すぎる場合はオーバーラップを調整
                excess = len(forward_seq) - 60
                overlap_length_adj = max(15, overlap_length - excess//2)
                left_overlap = template_seq[max(0, anneal_start - overlap_length_adj):anneal_start]
                right_overlap = template_seq[anneal_end:min(len(template_seq), anneal_end + overlap_length_adj)]
                forward_seq = left_overlap + annealing_region + right_overlap
                reverse_seq = self.reverse_complement(forward_seq)
            
            primers.append({
                'primer_number': primer_num,
                'codon_type': codon_type,
                'forward_name': f"{gene_name}_Sat_Fwd_P{primer_num}_{codon_type}",
                'forward_sequence': forward_seq,
                'forward_tm': round(self.calculate_tm(annealing_region), 1),
                'reverse_name': f"{gene_name}_Sat_Rev_P{primer_num}_{codon_type}",
                'reverse_sequence': reverse_seq,
                'reverse_tm': round(self.calculate_tm(annealing_region), 1),
                'combinations': combinations,
                'mixing_ratio': combinations,
                'primer_length': len(forward_seq),
                'annealing_length': len(annealing_region),
                'overlap_length': overlap_length
            })
        
        return primers
    
    def _design_double_nebuilder(self, template_seq: str, positions: List[int], gene_name: str, overlap_length: int) -> List[Dict]:
        """2箇所変異用NEBuilderプライマー"""
        primers = []
        pos1, pos2 = sorted(positions)
        
        # 両方の変異位置を含むアニーリング領域を設計
        center = (pos1 + pos2 + 3) // 2
        anneal_start, anneal_end = self._find_optimal_annealing_region(template_seq, center, min_length=24)
        
        # 変異位置がアニーリング領域に含まれるように調整
        anneal_start = min(anneal_start, pos1 - 3)
        anneal_end = max(anneal_end, pos2 + 6)
        
        primer_designs = [
            ('NDT', 'NDT', 1, 144), ('VHG', 'VHG', 2, 81), ('NDT', 'VHG', 3, 108),
            ('VHG', 'NDT', 4, 108), ('NDT', 'TGG', 5, 12), ('TGG', 'NDT', 6, 12),
            ('VHG', 'TGG', 7, 9), ('TGG', 'VHG', 8, 9), ('TGG', 'TGG', 9, 1)
        ]
        
        for codon1, codon2, primer_num, combinations in primer_designs:
            # アニーリング領域の構築
            left_part = template_seq[anneal_start:pos1]
            middle_part = template_seq[pos1 + 3:pos2]
            right_part = template_seq[pos2 + 3:anneal_end]
            annealing_region = left_part + codon1 + middle_part + codon2 + right_part
            
            # オーバーラップ領域
            left_overlap = template_seq[max(0, anneal_start - overlap_length):anneal_start]
            right_overlap = template_seq[anneal_end:min(len(template_seq), anneal_end + overlap_length)]
            
            forward_seq = left_overlap + annealing_region + right_overlap
            reverse_seq = self.reverse_complement(forward_seq)
            
            # 長さ調整
            if len(forward_seq) > 60:
                excess = len(forward_seq) - 60
                overlap_length_adj = max(15, overlap_length - excess//2)
                left_overlap = template_seq[max(0, anneal_start - overlap_length_adj):anneal_start]
                right_overlap = template_seq[anneal_end:min(len(template_seq), anneal_end + overlap_length_adj)]
                forward_seq = left_overlap + annealing_region + right_overlap
                reverse_seq = self.reverse_complement(forward_seq)
            
            primers.append({
                'primer_number': primer_num,
                'codon_types': [codon1, codon2],
                'forward_name': f"{gene_name}_Sat_Fwd_P{primer_num}_{codon1}{codon2}",
                'forward_sequence': forward_seq,
                'forward_tm': round(self.calculate_tm(annealing_region), 1),
                'reverse_name': f"{gene_name}_Sat_Rev_P{primer_num}_{codon1}{codon2}",
                'reverse_sequence': reverse_seq,
                'reverse_tm': round(self.calculate_tm(annealing_region), 1),
                'combinations': combinations,
                'mixing_ratio': combinations,
                'primer_length': len(forward_seq),
                'annealing_length': len(annealing_region),
                'overlap_length': overlap_length
            })
        
        return primers
    
    def calculate_statistics(self, sat_primers: List[Dict], vector_oligos: List[Dict], num_positions: int) -> Dict:
        """ライブラリー統計の計算（ベクターオリゴ含む）"""
        total_combinations = sum(p['combinations'] for p in sat_primers)
        unique_amino_acids = 20 ** num_positions
        coverage_95 = math.ceil(-unique_amino_acids * math.log(1 - 0.95))
        
        return {
            'total_sat_primers': len(sat_primers),
            'total_vector_oligos': len(vector_oligos),
            'total_combinations': total_combinations,
            'unique_amino_acids': unique_amino_acids,
            'coverage_95_percent': coverage_95,
            'library_efficiency': unique_amino_acids / total_combinations,
            'avg_primer_length': sum(p['primer_length'] for p in sat_primers) / len(sat_primers),
            'avg_annealing_tm': sum(p['forward_tm'] for p in sat_primers) / len(sat_primers)
        }

# Streamlit GUI
def main():
    st.set_page_config(
        page_title="22c-trick Saturation PCR Designer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    st.title("🧬 22c-trick Saturation PCR Designer")
    st.markdown("**NEBuilder最適化プライマー設計ツール**")
    
    # サイドバー：設定
    with st.sidebar:
        st.header("⚙️ 設定")
        
        gene_name = st.text_input("遺伝子名", value="MyGene", help="ファイル名とプライマー名に使用")
        
        st.markdown("### ベクター準備方法")
        vector_method = st.radio(
            "ベクターの準備方法を選択",
            ["PCR", "制限酵素"],
            help="PCR: ベクターをPCRで増幅\n制限酵素: ベクターを制限酵素で線状化"
        )
        
        st.markdown("### プライマー設計条件")
        overlap_length = st.slider("オーバーラップ長 (bp)", 15, 25, 20, help="NEBuilder用オーバーラップ領域")
        
        st.info(f"""
        **設計条件:**
        - アニーリング部: Tm ≥ 55°C, ≥ 18bp
        - オーバーラップ: {overlap_length}bp
        - 全長: 40bp推奨, 最大60bp
        - ベクター準備: {vector_method}
        
        **入力形式:**
        - DNA配列: 大文字・小文字・混在OK
        - 空白・改行は自動除去
        - 有効文字: A, T, G, C のみ
        """)
    
    # メイン画面：配列入力
    st.header("📝 配列入力")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        # サンプル配列ボタン
        col_sample1, col_sample2 = st.columns(2)
        
        with col_sample1:
            if st.button("🧪 EGFP サンプル配列"):
                sample_seq = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGAACTATGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
                st.session_state.template_seq_input = sample_seq
                st.session_state.template_seq = sample_seq
                st.rerun()
        
        with col_sample2:
            if st.button("🔤 混在形式テスト"):
                # 小文字・大文字・改行・スペースが混在したサンプル
                mixed_sample = """atggtgagcAAGGGCgaggagCTGTTCacc ggggtggtg
                cccatcCTGGTCGAGctggacggcGACGTAAAC
                ggccacaagTTCAGCgtgtccGGCGAGggc
                gagggcgatgccACCTACggcaagCTGACCctg
                aagttcatctgcACCACCggcaagCTGCCCgtg
                ccctggcccACCCTCgtgaccACCCTGaac
                tatggcgtgCAGTGCttcagcCGCTACccc
                gaccacatgAAGCAGcacgacTTCTTCaag
                tccgccatgCCCGAAggctacGTCCAGgag
                cgcaccatcTTCTTCaagGACGACggcaac
                tacaagaccCGCGCCgagGTGAAGttcgag
                ggcgacaccCTGGTGaaccGCATCgagctg
                aagggcatcGACTTCaagGAGGACggcaac
                atcctggggCACAAGctgGAGTACaactac
                aacagccacAACGTCtatatcATGGCCgac
                aagcagaagAACGGCatcAAGGTGaacttc
                aagatccgcCACAAcatcGAGGACggcagc
                gtgcagctcGCCGACcactACCAGcagaac
                acccccatcGGCGACggcCCCGTGctgctg
                cccgacaacCACTACctgAGCACCcagtcc
                gccctgagcAAAGACcccAACGAGaagcgc
                gatcacatgGTCCTGctgGAGTTCgtgacc
                gccgccgggATCACTctcGGCATGgacgag
                ctgtacaagtaa"""
                st.session_state.template_seq_input = mixed_sample
                st.session_state.template_seq = mixed_sample.upper().replace(' ', '').replace('\n', '').replace('\r', '')
                st.rerun()
        
        template_seq_input = st.text_area(
            "テンプレートDNA配列 (CDS全長)",
            value=st.session_state.get('template_seq_input', ''),
            height=150,
            help="ATGから終止コドンまでの全CDS配列を入力（大文字・小文字・混在OK）"
        )
        
        if template_seq_input:
            # DNA配列の前処理
            template_seq = template_seq_input.upper().replace(' ', '').replace('\n', '').replace('\r', '')
            st.session_state.template_seq_input = template_seq_input
            st.session_state.template_seq = template_seq
            
            # 無効文字チェック
            valid_bases = set('ATGC')
            invalid_chars = set(template_seq) - valid_bases
            if invalid_chars:
                st.error(f"⚠️ 無効な文字が含まれています: {', '.join(sorted(invalid_chars))}")
                st.info("有効な文字: A, T, G, C (大文字・小文字・混在OK)")
            else:
                if template_seq_input != template_seq:
                    st.info(f"✅ 自動変換: 大文字に統一、空白・改行を除去")
    
    with col2:
        if template_seq_input:
            # 前処理済み配列を取得
            template_seq = st.session_state.get('template_seq', '')
            
            if template_seq:
                seq_length = len(template_seq)
                aa_length = seq_length // 3
                st.metric("配列長", f"{seq_length} bp")
                st.metric("アミノ酸数", f"{aa_length} aa")
                
                # 配列の妥当性チェック
                valid_bases = set('ATGC')
                invalid_chars = set(template_seq) - valid_bases
                
                if invalid_chars:
                    st.error("❌ 無効文字あり")
                elif seq_length % 3 != 0:
                    st.warning("⚠️ 3の倍数でない")
                elif not template_seq.startswith('ATG'):
                    st.warning("⚠️ ATGで始まらない")
                elif template_seq[:-3].find('TAA') != -1 or template_seq[:-3].find('TAG') != -1 or template_seq[:-3].find('TGA') != -1:
                    st.warning("⚠️ 途中にSTOPコドン")
                else:
                    st.success("✅ 配列形式OK")
    
    # アミノ酸位置選択
    template_seq = st.session_state.get('template_seq', '')
    if template_seq and len(template_seq) % 3 == 0 and template_seq.startswith('ATG'):
        # 無効文字チェック
        valid_bases = set('ATGC')
        invalid_chars = set(template_seq) - valid_bases
        
        if not invalid_chars:  # 有効な配列の場合のみ処理
            st.header("🎯 変異位置選択")
            
            designer = OptimizedSaturationDesigner()
            protein_seq = designer.translate_dna(template_seq)
            aa_length = len(protein_seq)
            
            # アミノ酸配列表示（10個ずつ）
            st.markdown("### アミノ酸配列")
            for i in range(0, len(protein_seq), 50):
                line_start = i + 1
                line_end = min(i + 50, len(protein_seq))
                line_seq = protein_seq[i:line_end]
                
                # 位置番号を表示
                pos_str = "".join([f"{j:>2}" if (j-1) % 10 == 0 else "  " for j in range(line_start, line_end + 1)])
                seq_str = "".join([f"{aa:>2}" for aa in line_seq])
                
                st.text(f"{line_start:>3}-{line_end:>3}: {seq_str}")
        
        # アミノ酸配列表示（10個ずつ）
        st.markdown("### アミノ酸配列")
        for i in range(0, len(protein_seq), 50):
            line_start = i + 1
            line_end = min(i + 50, len(protein_seq))
            line_seq = protein_seq[i:line_end]
            
            # 位置番号を表示
            pos_str = "".join([f"{j:>2}" if (j-1) % 10 == 0 else "  " for j in range(line_start, line_end + 1)])
            seq_str = "".join([f"{aa:>2}" for aa in line_seq])
            
            st.text(f"{line_start:>3}-{line_end:>3}: {seq_str}")
        
        # アミノ酸位置選択
        st.subheader("🎯 変異位置選択")
        
        col1, col2 = st.columns(2)
        
        with col1:
            position_1 = st.selectbox(
                "1番目の変異位置",
                options=list(range(1, aa_length + 1)),
                format_func=lambda x: f"{x}: {protein_seq[x-1]}",
                help="変異を導入したいアミノ酸位置"
            )
        
        with col2:
            enable_second = st.checkbox("2箇所目も変異させる")
            if enable_second:
                position_2 = st.selectbox(
                    "2番目の変異位置",
                    options=[i for i in range(1, aa_length + 1) if i != position_1],
                    format_func=lambda x: f"{x}: {protein_seq[x-1]}",
                    help="2番目の変異位置"
                )
            else:
                position_2 = None
        
        # 制限酵素サイト選択（制限酵素モードの場合）
        selected_enzyme = None
        selected_position = None
        
        if vector_method == "制限酵素":
            st.subheader("🔪 制限酵素サイト選択")
            
            restriction_sites = designer.find_restriction_sites(template_seq)
            
            if restriction_sites:
                # 制限酵素サイト表示
                st.markdown("### 検出された制限酵素サイト")
                
                site_options = []
                for enzyme, positions in restriction_sites.items():
                    for pos in positions:
                        site_options.append({
                            'enzyme': enzyme,
                            'position': pos,
                            'display': f"{enzyme} at {pos}bp ({designer.restriction_enzymes[enzyme]['sequence']})"
                        })
                
                if site_options:
                    selected_site = st.selectbox(
                        "使用する制限酵素サイトを選択",
                        options=range(len(site_options)),
                        format_func=lambda x: site_options[x]['display'],
                        help="ベクター線状化に使用する制限酵素サイト"
                    )
                    
                    selected_enzyme = site_options[selected_site]['enzyme']
                    selected_position = site_options[selected_site]['position']
                    
                    st.success(f"✅ 選択: {selected_enzyme} at {selected_position}bp")
                else:
                    st.warning("⚠️ 制限酵素サイトが見つかりませんでした")
            else:
                st.warning("⚠️ 認識可能な制限酵素サイトがありません")
                st.info("💡 PCRモードを使用することをお勧めします")
        
        # 設計実行
        design_ready = True
        if vector_method == "制限酵素" and (selected_enzyme is None or selected_position is None):
            design_ready = False
        
        if st.button("🚀 プライマー設計実行", type="primary", disabled=not design_ready):
            mutation_positions = [position_1]
            if position_2:
                mutation_positions.append(position_2)
            
            with st.spinner("設計中..."):
                # DNA位置に変換（0-based）
                dna_positions = [(pos - 1) * 3 for pos in mutation_positions]
                
                # プライマー設計（新しいAPI）
                sat_primers, vector_oligos = designer.design_nebuilder_primers(
                    template_seq, dna_positions, gene_name, overlap_length,
                    vector_method, selected_enzyme, selected_position
                )
                
                # 統計計算
                stats = designer.calculate_statistics(sat_primers, vector_oligos, len(mutation_positions))
                
                # セッション状態に保存
                st.session_state.sat_primers = sat_primers
                st.session_state.vector_oligos = vector_oligos
                st.session_state.stats = stats
                st.session_state.mutation_positions = mutation_positions
                st.session_state.vector_method = vector_method
    
    # 結果表示
    if 'sat_primers' in st.session_state:
        st.header("📊 設計結果")
        
        # 統計情報
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("Saturation", f"{st.session_state.stats['total_sat_primers']} プライマー")
        with col2:
            st.metric("Vector", f"{st.session_state.stats['total_vector_oligos']} オリゴ")
        with col3:
            st.metric("ライブラリー", f"{st.session_state.stats['total_combinations']} 組合せ")
        with col4:
            st.metric("95%カバレッジ", f"{st.session_state.stats['coverage_95_percent']} クローン")
        with col5:
            st.metric("効率", f"{st.session_state.stats['library_efficiency']:.1%}")
        
        # ベクター準備方法の表示
        vector_method = st.session_state.get('vector_method', 'PCR')
        if vector_method == "PCR":
            st.info("🧪 **ベクター準備**: PCRによる増幅")
        else:
            st.info("🔪 **ベクター準備**: 制限酵素による線状化")
        
        # Saturation Mutagenesis プライマー詳細
        st.subheader("🧬 Saturation Mutagenesis プライマー")
        
        # Saturationプライマーテーブル作成
        sat_primer_data = []
        for primer in st.session_state.sat_primers:
            sat_primer_data.extend([
                {
                    'Name': primer['forward_name'],
                    'Sequence': primer['forward_sequence'],
                    'Length': primer['primer_length'],
                    'Annealing Tm': primer['forward_tm'],
                    'Type': 'Forward',
                    'Mixing Ratio': primer['mixing_ratio'],
                    'Primer #': primer['primer_number']
                },
                {
                    'Name': primer['reverse_name'],
                    'Sequence': primer['reverse_sequence'],
                    'Length': primer['primer_length'],
                    'Annealing Tm': primer['reverse_tm'],
                    'Type': 'Reverse',
                    'Mixing Ratio': primer['mixing_ratio'],
                    'Primer #': primer['primer_number']
                }
            ])
        
        sat_df = pd.DataFrame(sat_primer_data)
        st.dataframe(sat_df, use_container_width=True)
        
        # Vector用オリゴ詳細
        if st.session_state.vector_oligos:
            st.subheader(f"🔧 Vector用オリゴ ({vector_method})")
            
            vector_data = []
            for oligo in st.session_state.vector_oligos:
                vector_data.append({
                    'Name': oligo['name'],
                    'Sequence': oligo['sequence'],
                    'Length': oligo['length'],
                    'Tm': oligo['tm'],
                    'Type': oligo['type'],
                    'Purpose': oligo['purpose'],
                    'Description': oligo['description']
                })
            
            vector_df = pd.DataFrame(vector_data)
            st.dataframe(vector_df, use_container_width=True)
        
        # 混合比率表示
        st.subheader("⚗️ プライマー混合比率")
        forward_primers = [p for p in st.session_state.sat_primers]
        mix_info = []
        for primer in forward_primers:
            mix_info.append({
                'プライマー番号': primer['primer_number'],
                'Forward': primer['forward_name'],
                'Reverse': primer['reverse_name'],
                '混合比率': primer['mixing_ratio']
            })
        
        mix_df = pd.DataFrame(mix_info)
        st.dataframe(mix_df, use_container_width=True)
        
        # 実験プロトコル
        st.subheader("🔬 実験プロトコル")
        
        if vector_method == "PCR":
            st.markdown("""
            **PCRによるベクター準備:**
            1. Vector用PCRオリゴを使用してベクターを増幅
            2. PCR産物を精製
            3. Saturation primersでinsert増幅
            4. NEBuilder HiFi DNA Assemblyで組み立て
            """)
        else:
            st.markdown("""
            **制限酵素によるベクター準備:**
            1. 選択した制限酵素でベクターを線状化
            2. ベクターを精製（アルカリホスファターゼ処理推奨）
            3. Saturation primersでinsert増幅
            4. NEBuilder HiFi DNA Assemblyで組み立て
            """)
        
        # ファイルダウンロード
        st.subheader("💾 ファイルダウンロード")
        
        # 全オリゴの統合
        all_oligos = sat_primer_data.copy()
        if st.session_state.vector_oligos:
            for oligo in st.session_state.vector_oligos:
                all_oligos.append({
                    'Name': oligo['name'],
                    'Sequence': oligo['sequence'],
                    'Length': oligo['length'],
                    'Annealing Tm': oligo.get('tm', 'N/A'),
                    'Type': oligo['type'],
                    'Mixing Ratio': 'N/A',
                    'Primer #': 'Vector'
                })
        
        all_df = pd.DataFrame(all_oligos)
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv_data = all_df.to_csv(index=False)
            st.download_button(
                "📄 全オリゴテーブル (CSV)",
                csv_data,
                f"{gene_name}_{vector_method}_all_oligos.csv",
                "text/csv"
            )
        
        with col2:
            # 受注用フォーマット
            order_text = "Oligo Name\tSequence\n"
            for _, row in all_df.iterrows():
                if 'Sequence' in row and len(str(row['Sequence'])) > 10:  # 有効な配列のみ
                    order_text += f"{row['Name']}\t{row['Sequence']}\n"
            
            st.download_button(
                "📋 受注用フォーマット (TXT)",
                order_text,
                f"{gene_name}_{vector_method}_order.txt",
                "text/plain"
            )
        
        with col3:
            # 統計レポート
            stats_text = f"""22c-trick Saturation PCR Design Report
Gene: {gene_name}
Vector Method: {vector_method}
Mutation positions: {', '.join(map(str, st.session_state.mutation_positions))}

Statistics:
- Saturation primers: {st.session_state.stats['total_sat_primers']}
- Vector oligos: {st.session_state.stats['total_vector_oligos']}
- Library combinations: {st.session_state.stats['total_combinations']}
- Unique amino acids: {st.session_state.stats['unique_amino_acids']}
- 95% coverage: {st.session_state.stats['coverage_95_percent']} clones
- Library efficiency: {st.session_state.stats['library_efficiency']:.3f}
- Average primer length: {st.session_state.stats['avg_primer_length']:.1f} bp
- Average annealing Tm: {st.session_state.stats['avg_annealing_tm']:.1f}°C

Mixing Ratios (Saturation Primers):
"""
            for primer in forward_primers:
                stats_text += f"P{primer['primer_number']}: {primer['mixing_ratio']} parts\n"
            
            if vector_method == "制限酵素":
                stats_text += f"\nRestriction enzyme protocol:\n"
                stats_text += f"- Use selected restriction enzyme for vector linearization\n"
                stats_text += f"- Treat with alkaline phosphatase if needed\n"
            else:
                stats_text += f"\nPCR protocol:\n"
                stats_text += f"- Use Vector PCR oligos for vector amplification\n"
                stats_text += f"- Purify PCR products before assembly\n"
            
            st.download_button(
                "📈 実験レポート (TXT)",
                stats_text,
                f"{gene_name}_{vector_method}_report.txt",
                "text/plain"
            )

if __name__ == "__main__":
    main()
