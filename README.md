# Genome synteny and duplication

用于 Codex 的基因组共线性与复制分析 Skill。以完整代表蛋白组和基因坐标建立可复用的共线性证据，再叠加任意关注基因家族；保留原始基因 ID、参数、输入输出哈希和图形来源表。

## 功能与范围

- 基因组内部的严格串联复制候选、近邻复制候选和 WGD/segmental 共线性候选。
- 亚基因组、物种或不同组装之间的共线性区块、直接锚定基因对及块内所有基因。
- 新增基因家族时复用已有全基因组分析，不按家族重新选择代表蛋白。
- 使用 Circos 绘制环形图、JCVI 绘制双基因组染色体图，以及 Matplotlib 绘制点图和局部基因顺序图；附实际绘制的基因、连线、图注和省略数量。

Skill 位于 [`skills/genome-synteny-duplication`](skills/genome-synteny-duplication/SKILL.md)。共线性、邻接和家族归属分别保留；这些证据不能单独证明酶功能、特定 WGD 事件或基因丢失。最终发表尺寸下仍须检查图形。

## 安装到 Codex

克隆仓库后，在仓库根目录操作：

```bash
git clone https://github.com/yz3394/Genome-synteny-duplication.git
cd Genome-synteny-duplication
```

将 `skills/genome-synteny-duplication` 整个目录复制到 Codex 的 `skills` 目录。**已有同名本地 Skill 时，先保存并比较本地修改，不要直接覆盖。** 下列命令遇到已有目录会跳过复制：

```bash
skill_target="${CODEX_HOME:-$HOME/.codex}/skills/genome-synteny-duplication"
if [ -e "$skill_target" ]; then
  echo "目标已存在：请先保存并比较本地修改。未执行复制。"
else
  mkdir -p "${CODEX_HOME:-$HOME/.codex}/skills"
  cp -R skills/genome-synteny-duplication "$skill_target"
fi
```

在能读取该 Skill 的 Codex 会话中使用 `$genome-synteny-duplication`，提供匹配的 GFF3、蛋白 FASTA、基因组长度或组装 FASTA，以及已确认的家族成员列表。

## 运行环境

| 用途 | 依赖 |
|---|---|
| 输入准备、证据整理、JSON 配置 | Python ≥3.10，标准库 |
| 全蛋白组搜索 | DIAMOND，或 BLAST+ 的 `blastp` 与 `makeblastdb` |
| 共线性推断 | MCScanX；原生基因分类另需 `duplicate_gene_classifier` |
| 局部图、点图 | Matplotlib |
| 双基因组染色体图 | JCVI 及其运行环境 |
| 环形图 | Circos、Matplotlib、Pillow |
| YAML 配置（可选） | PyYAML |

本仓库不包含这些外部软件、数据库或基因组。开发验证在 macOS 完成；其他平台和软件版本须检查实际输出参数及相关测试。开发使用的 MCScanX 提交为 `0956cb8f900c152e2be8ba3196829e50ab454b94`；新的环境应记录自身版本。

## 基本用法

在仓库根目录执行。先复制并编辑 [`analysis_config.example.json`](skills/genome-synteny-duplication/assets/analysis_config.example.json)，替换基因组信息、输入位置和软件路径；输入路径相对于配置文件解析。

```bash
python3 skills/genome-synteny-duplication/scripts/synteny_workflow.py prepare --config /path/to/analysis.json
python3 skills/genome-synteny-duplication/scripts/synteny_workflow.py build --config /path/to/analysis.json
python3 skills/genome-synteny-duplication/scripts/synteny_workflow.py overlay --run /path/to/run_manifest.json --families /path/to/families.tsv
```

默认起始方案为 DIAMOND sensitive、蛋白 E≤1e−5、报告上限 50；MCScanX 使用每个查询的前 5 个不同非自身目标。串联候选保留较宽的同源命中池。参数是起点，应根据命中饱和、序列差异和重点位点进行必要的敏感性检查。

- [输入格式与证据表](skills/genome-synteny-duplication/references/inputs-and-evidence.md)
- [方法、参数和原始文献](skills/genome-synteny-duplication/references/methods-and-parameters.md)
- [绘图命令、配置与检查](skills/genome-synteny-duplication/references/figures-and-qa.md)

## 验证范围

2026-09-28 的开发验证包括咖啡 CC/CE 与长春花 v3 的五组全蛋白组比较、旧结果回归、参数与代表转录本敏感性检查、六组区域性搜索配置，以及 Circos/JCVI/局部图的实际渲染。详情见[验证范围](skills/genome-synteny-duplication/references/coffee-regression.md)；项目原始数据和完整运行结果未随仓库发布。

当时完整启用环境通过 **52 项测试**。在新的环境可运行：

```bash
python3 -m unittest discover -s skills/genome-synteny-duplication/tests
SYNTENY_MCSCANX=/absolute/path/to/MCScanX python3 -m unittest discover -s skills/genome-synteny-duplication/tests
```

未设置 `SYNTENY_MCSCANX` 会跳过 3 项真实引擎测试；未安装 Matplotlib 另跳过 2 项渲染测试。历史验证不代表新的数据、平台或依赖版本已经验证，也不代表区域性 BLASTP 检查已证明全蛋白组搜索等价。

## 本地更新与 GitHub 同步

仓库所有者已授权将本地安装的 Skill 更新同步至此仓库。同步方向、每小时检查的运行条件、验证门槛及冲突处理见 [MAINTENANCE.md](MAINTENANCE.md)。文件复制、Git 提交和 GitHub 推送是独立步骤；同步脚本自身不执行推送。
