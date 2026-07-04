# Changelog

## 1.1.0 (2026-07-05)
- 修复 star_read_command 参数命名与 Snakemake 兼容性
- STAR 规则添加 `--limitBAMsortRAM` 内存资源声明
- samtools index 添加多线程支持 (`-@`)
- MultiQC 精确扫描指定目录 (results/fastqc, results/trimmed)
- DESeq2 输出 sessionInfo() 确保可复现性
- 添加 test data 下载脚本 (scripts/download_test_data.sh)
- 添加 GitHub Actions CI 工作流 (.github/workflows/test.yml)
- 添加 Dockerfile 与 .dockerignore 容器支持
- 更新 reference/README.md 下载和索引构建说明
- 更新 .gitignore 排除参考基因组子目录
- 移除空 report/ 目录，data/README.md 中英双语
- README 添加 CI badge、Docker 使用说明、测试数据步骤

## 1.0.0 (2026-07-04)
- 完整 RNA-seq 分析流程
- 输入校验、FastQC 质控、STAR 比对、featureCounts 定量
- 表达矩阵合并、DESeq2 差异分析
- 样本表驱动配置，支持双端/单端
- 可选 fastp 修剪
- Conda 分步环境管理
- 输入校验与错误提示
- chr22 测试参考基因组
