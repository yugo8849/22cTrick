# 🧬 22c-trick Saturation PCR Designer

**NEBuilder最適化プライマー設計ツール**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://saturation-pcr-designer.streamlit.app)

## 🎯 特徴

- **22c-trick手法**による高効率saturation mutagenesis
- **STOPコドン完全回避**で品質保証
- **NEBuilder最適化**プライマー設計
- **視覚的GUI**で直感的操作
- **1-2箇所の同時変異**対応
- **大文字・小文字・混在**入力対応

## 🚀 オンライン使用

[**🌐 アプリを開く**](https://saturation-pcr-designer.streamlit.app)

## 📊 プライマー設計条件

| 項目 | 条件 |
|------|------|
| アニーリング部分 | Tm ≥ 55°C, ≥ 18bp |
| オーバーラップ | 20bp固定 |
| 全長 | 40bp推奨, 最大60bp |
| 変異方式 | 22c-trick (NDT/VHG/TGG) |

## 📖 使用方法

### 1. **配列入力**
```
# どの形式でもOK
ATGGTGAGCAAG...     # 大文字
atggtgagcaag...     # 小文字  
ATGgtgAGCaag...     # 混在
ATG GTG AGC AAG...  # スペース・改行入り
```

### 2. **変異位置選択**
- アミノ酸配列表示から位置を選択
- 1箇所または2箇所の変異をサポート

### 3. **プライマー設計**
- ワンクリックで全プライマーを設計
- 統計情報とライブラリー効率を表示

### 4. **結果ダウンロード**
- **詳細テーブル** (.csv): 全プライマー情報
- **受注用フォーマット** (.txt): オリゴ合成業者用
- **統計レポート** (.txt): 実験計画用

## 🔬 実験手順

1. **オリゴ注文**: 受注用ファイルを合成業者にコピペ
2. **プライマー混合**: 表示された比率で混合
3. **PCR**: QuikChange風プロトコル（アニーリング60°C推奨）
4. **DpnI消化**: テンプレート除去
5. **形質転換**: 統計レポートのクローン数を参考
6. **スクリーニング**: 22c-trick効率で工数削減

## 📈 効率比較

| 手法 | コドン数 | 2箇所変異 | STOPコドン | 効率 |
|------|----------|-----------|------------|------|
| NNN | 64 | 4,096 | あり(4.7%) | 低 |
| NNK/NNS | 32 | 1,024 | あり(3.1%) | 中 |
| **22c-trick** | **22** | **484** | **なし** | **高** |

## 🛠️ ローカル実行

```bash
# Conda環境作成
conda create -n saturation_pcr python=3.9 -y
conda activate saturation_pcr

# 依存関係インストール
pip install -r requirements.txt

# アプリ実行
streamlit run app.py
```

## 📚 科学的根拠

このツールは以下の査読済み論文に基づいています：

> Kille, S., et al. "Reducing codon redundancy and screening effort of combinatorial protein libraries created by saturation mutagenesis." *ACS Synthetic Biology* (2013)

**主な利点:**
- ライブラリーサイズ50%削減
- スクリーニング工数大幅短縮
- STOPコドン完全排除
- 統計的最適化

## 🎮 サンプル機能

- **🧪 EGFP サンプル**: 標準的なGFP配列でのテスト
- **🔤 混在形式テスト**: 大文字・小文字・改行混在のテスト

## 🔧 技術仕様

- **フレームワーク**: Streamlit
- **言語**: Python 3.9+
- **依存関係**: pandas, numpy
- **デプロイ**: Streamlit Cloud
- **ライセンス**: MIT

## 📊 プロジェクト統計

- コード行数: ~500行
- 対応生物: すべて（universal genetic code）
- プライマー長: 25-60bp
- 変異位置: 1-2箇所同時

## 🤝 貢献

Issues、Pull Requests歓迎です！

### 開発者向け

```bash
# リポジトリフォーク後
git clone https://github.com/[YOUR_USERNAME]/saturation-pcr-designer.git
cd saturation-pcr-designer

# 開発環境セットアップ
conda env create -f environment.yml
conda activate saturation_pcr

# ブランチ作成
git checkout -b feature/new-feature

# 変更後
git commit -m "Add: new feature"
git push origin feature/new-feature
```

## 📞 サポート

- **Issues**: バグ報告・機能要望
- **Discussions**: 使用方法・質問
- **Email**: 研究協力のお問い合わせ

## 📝 ライセンス

MIT License - 研究・商用利用自由

## 🙏 謝辞

- 22c-trick手法の開発者の皆様
- Streamlitコミュニティ
- 分子生物学研究者の皆様

---

**🧬 高効率なsaturation mutagenesis実験をお楽しみください！**