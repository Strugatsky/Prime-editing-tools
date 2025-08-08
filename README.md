# Prime-editing-tools
A collection of tools intended to automate the prime editing workflow of the GEML at ETH. These tools were written for my specific workflow and not in a robust, generalizable way. However, it should be fairly simple to make them work with most kinds of data structures.
## Original Workflow
![Original Workflow](figures/original_workflow.drawio.pdf)
## New Workflow
![New Workflow](figures/new_workflow.drawio.pdf)
## Tools
### Crispresso run generator
```
python crispresso_run_generator.py  [-h]                            \
                                    --seq_folder SEQ_FOLDER         \
                                    --database DATABASE             \
                                    --amplicon_file AMPLICON_FILE   \
                                    --scaffold_file SCAFFOLD_FILE   \
                                    [--output OUTPUT]
```
### Oligo order generator
```
python oligo_order_generator.py [-h]                \
                                [--prefix PREFIX]   \
                                [--use_scaffold_2]  \
                                input_file          \
                                output_file         \
                                db_file
```
### Quantification analysis
```
python quantification_analysis.py   [-h]        \
                                    tsv_file    \
                                    db_file
```
