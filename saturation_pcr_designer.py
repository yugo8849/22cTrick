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
    
    def translate_dna(self, dna_seq: str) -> str:
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
                                gene_name: str = "Gene", overlap_length: int = 20) -> List[Dict]:
        """NEBuilder最適化プライマー設計"""
        
        if len(mutation_positions) == 1:
            return self._design_single_nebuilder(template_seq, mutation_positions[0], gene_name, overlap_length)
        elif len(mutation_positions) == 2:
            return self._design_double_nebuilder(template_seq, mutation_positions, gene_name, overlap_length)
        else:
            raise ValueError("1〜2箇所の変異位置のみ対応")
    
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
    
    def calculate_statistics(self, primers: List[Dict], num_positions: int) -> Dict:
        """統計計算"""
        total_combinations = sum(p['combinations'] for p in primers)
        unique_amino_acids = 20 ** num_positions
        coverage_95 = math.ceil(-unique_amino_acids * math.log(1 - 0.95))
        
        return {
            'total_primers': len(primers),
            'total_combinations': total_combinations,
            'unique_amino_acids': unique_amino_acids,
            'coverage_95_percent': coverage_95,
            'library_efficiency': unique_amino_acids / total_combinations,
            'avg_primer_length': sum(p['primer_length'] for p in primers) / len(primers),
            'avg_annealing_tm': sum(p['forward_tm'] for p in primers) / len(primers)
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
        
        st.markdown("### プライマー設計条件")
        overlap_length = st.slider("オーバーラップ長 (bp)", 15, 25, 20, help="NEBuilder用オーバーラップ領域")
        
        st.info("""
        **設計条件:**
        - アニーリング部: Tm ≥ 55°C, ≥ 18bp
        - オーバーラップ: 20bp固定
        - 全長: 40bp推奨, 最大60bp
        
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
        
        st.markdown("### 変異位置選択")
        
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
        
        # 設計実行
        if st.button("🚀 プライマー設計実行", type="primary"):
            mutation_positions = [position_1]
            if position_2:
                mutation_positions.append(position_2)
            
            with st.spinner("設計中..."):
                # DNA位置に変換（0-based）
                dna_positions = [(pos - 1) * 3 for pos in mutation_positions]
                
                # プライマー設計
                primers = designer.design_nebuilder_primers(
                    template_seq, dna_positions, gene_name, overlap_length
                )
                
                # 統計計算
                stats = designer.calculate_statistics(primers, len(mutation_positions))
                
                # セッション状態に保存
                st.session_state.primers = primers
                st.session_state.stats = stats
                st.session_state.mutation_positions = mutation_positions
    
    # 結果表示
    if 'primers' in st.session_state:
        st.header("📊 設計結果")
        
        # 統計情報
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("プライマー数", st.session_state.stats['total_primers'])
        with col2:
            st.metric("ライブラリー組合せ", st.session_state.stats['total_combinations'])
        with col3:
            st.metric("95%カバレッジ", f"{st.session_state.stats['coverage_95_percent']} クローン")
        with col4:
            st.metric("ライブラリー効率", f"{st.session_state.stats['library_efficiency']:.1%}")
        
        # プライマー詳細
        st.subheader("🧪 プライマー詳細")
        
        # テーブル作成
        primer_data = []
        for primer in st.session_state.primers:
            primer_data.extend([
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
        
        df = pd.DataFrame(primer_data)
        st.dataframe(df, use_container_width=True)
        
        # 混合比率表示
        st.subheader("⚗️ プライマー混合比率")
        forward_primers = [p for p in st.session_state.primers]
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
        
        # ファイルダウンロード
        st.subheader("💾 ファイルダウンロード")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv_data = df.to_csv(index=False)
            st.download_button(
                "📄 詳細テーブル (CSV)",
                csv_data,
                f"{gene_name}_primers.csv",
                "text/csv"
            )
        
        with col2:
            # 受注用フォーマット
            order_text = "Oligo Name\tSequence\n"
            for _, row in df.iterrows():
                order_text += f"{row['Name']}\t{row['Sequence']}\n"
            
            st.download_button(
                "📋 受注用フォーマット (TXT)",
                order_text,
                f"{gene_name}_order.txt",
                "text/plain"
            )
        
        with col3:
            # 統計レポート
            stats_text = f"""22c-trick Saturation PCR Design Report
Gene: {gene_name}
Mutation positions: {', '.join(map(str, st.session_state.mutation_positions))}

Statistics:
- Total primers: {st.session_state.stats['total_primers']}
- Library combinations: {st.session_state.stats['total_combinations']}
- Unique amino acids: {st.session_state.stats['unique_amino_acids']}
- 95% coverage: {st.session_state.stats['coverage_95_percent']} clones
- Library efficiency: {st.session_state.stats['library_efficiency']:.3f}
- Average primer length: {st.session_state.stats['avg_primer_length']:.1f} bp
- Average annealing Tm: {st.session_state.stats['avg_annealing_tm']:.1f}°C

Mixing Ratios:
"""
            for primer in forward_primers:
                stats_text += f"P{primer['primer_number']}: {primer['mixing_ratio']} parts\n"
            
            st.download_button(
                "📈 統計レポート (TXT)",
                stats_text,
                f"{gene_name}_stats.txt",
                "text/plain"
            )

if __name__ == "__main__":
    main()
