# 结果如何回到Mac

## 回传内容

所有结果放本交接目录的`results/`，图放`01_图形`、表放`02_表格`、对象放`03_对象`、日志/检查放`04_日志和检查`，新增后处理代码放`code/`。

必须有：

- `RUN_SUMMARY.md`：按TASK_SPEC逐项标明完成/失败/未做，实际细胞范围、计算耗时、软件版本、警告和未解决问题。
- `01_逐图阅读指南.md`：相对路径嵌入每张PNG，链接PDF及源表；逐图写实际发现，不只解释读图。
- 原始DDRTree断点、Monocle2主/替代根RDS、Slingshot对象与全部逐细胞/路径/树表。
- 输入/代码/参数/环境绑定，数值检查、全部PNG/PDF视觉记录、stdout/stderr与退出状态。
- 所有后处理代码、完整Python依赖清单、R sessionInfo；任何修复的变更说明。

不要求把1GB原始输入重复回传。不要上传系统Keychain、.env、SSH密钥、访问令牌、整机目录或版权PDF。

## 代码和大文件分开

1. 在本仓库新建明确命名的回传分支，例如`pc-results-20261022`，先核对当前分支/工作树，不覆盖他人修改。
2. 提交必要代码、说明和小型摘要。`data/work/results/return`默认在gitignore里，大对象不要强行git add。
3. 使用入口打包，例如：

```bash
python Monocle2_PC.py pack-return --label pc_results_20261022
```

上面的日期只是命名示例，实际执行时使用当天或新版本标签，不覆盖已有目录。

4. 将`return/<标签>/RETURN_MANIFEST.json`、`RESULT_MANIFEST.json`、RUN_SUMMARY和清单列出的结果附件上传到一个**新的结果Release**，并把Release指向对应代码提交。每个附件小于2GiB。入口超过1.5GiB时只做传输文件分片，不拆开细胞拟合。
5. 发布后读取远端附件列表，核对文件名、字节数和服务端SHA256。README/运行摘要写明结果Release链接、代码commit、输入Release版本和是否通过验收。没有权限时先保存本地完整包，再向用户索取协作权限；不能删除结果或使用其他账号秘密。

示例命令形态（先替换经核实的标签、提交和附件列表）：

```text
gh release create <新的结果标签> <清单列出的附件> --repo COSMOS340/exchange --target <结果代码提交SHA> --title <结果版本名称> --notes-file results/RUN_SUMMARY.md
```

不要上传巨大的本地`results.zip`和它的分片两份。无分片时传zip；有分片时只传清单列出的part文件及两份manifest。不要使用`--clobber`覆盖输入Release或旧结果。

## Mac接收验收

Mac从确切的结果Release下载，先核对RETURN_MANIFEST的每个附件SHA和字节数；若分片，按清单次序二进制拼接，复核整包archive_sha256，再解压到新的接收目录。按RESULT_MANIFEST逐文件核验。

随后检查150,917细胞、CNV冻结值、有效根/树/路径、动态统计单位和全部图/PDF。仅在数值与视觉验收通过后更新第14步正式入口；之前的Mac输入/旧结果不删除。PC端“完成”不能代替Mac复核。

回传结束在聊天提供仓库分支、commit、Release链接和简短完成/失败清单。GitHub不会自动通知本地GPT取结果，用户需告诉Mac端“PC已回传”，再进行下载与复核。

[GitHub Release附件限制](https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases)
