import sys

def append_issue(file_path, priority, issue_line):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    header = f"## {priority}"
    target_idx = -1
    for i, line in enumerate(lines):
        if line.strip() == header:
            target_idx = i
            break
    
    if target_idx != -1:
        lines.insert(target_idx + 1, issue_line + "\n")
    else:
        lines.append("\n" + header + "\n" + issue_line + "\n")
        
    with open(file_path, 'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    file_path = "ISSUES.md"
    issues = [
        ("P1 Critical", "- [ ] #77 [P1] CanSNPer2/CanSNPer2/modules/NewickTree.py:288 — UnboundLocalError in draw_ete3_tree if ete3.Tree fails (sweep)"),
        ("P2 High", "- [ ] #78 [P2] CanSNPer2/CanSNPer2/modules/NewickTree.py:387 — IndexError if tree.search_nodes returns no results in draw_ete3_tree (sweep)"),
        ("P2 High", "- [ ] #79 [P2] CanSNPer2/CanSNPer2/modules/ModifyTree.py:33 — Integer division in progressBar causes progress bar to stay at 0% (sweep)"),
        ("P2 High", "- [ ] #80 [P2] CanSNPer2/CanSNPer2/modules/ModifyTree.py:477 — update_genomes might be setting set_value to the literal string \"id\" (sweep)"),
        ("P2 High", "- [ ] #81 [P2] CanSNPer2/CanSNPer2/modules/ModifyTree.py:267 — IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        ("P2 High", "- [ ] #82 [P2] CanSNPer2/CanSNPer2/modules/ModifyTree.py:284 — IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        ("P2 High", "- [ ] #83 [P2] CanSNPer2/CanSNPer2/modules/ModifyTree.py:368 — IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        ("P2 High", "- [ ] #84 [P2] CanSNPer2/CanSNPer2/modules/ModifyTree.py:370 — IndexError/ValueError in file_mod if modification file is malformed (sweep)"),
        ("P2 High", "- [ ] #85 [P2] CanSNPer2/CanSNPer2/modules/runMauve.py:81 — fragile command parsing with split(\" \") in run_mauve (sweep)"),
        ("P2 High", "- [ ] #86 [P2] CanSNPer2/CanSNPer2/modules/runMauve.py:142 — cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
        ("P2 High", "- [ ] #87 [P2] CanSNPer2/CanSNPer2/modules/runMauve.py:143 — cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
        ("P3 Medium", "- [ ] #88 [P3] CanSNPer2/CanSNPer2/modules/ModifyDatabase.py:77 — IndexError in load_annotation_file if data row is shorter than headers (sweep)"),
        ("P3 Medium", "- [ ] #89 [P3] CanSNPer2/CanSNPer2/modules/ModifyDatabase.py:107 — IndexError in load_genome_reference_file if data row is shorter than headers (sweep)"),
        ("P3 Medium", "- [ ] #90 [P3] CanSNPer2/CanSNPer2/modules/ModifyTree.py:292 — KeyError in update_annotations if name is not in snplist (sweep)"),
        ("P3 Medium", "- [ ] #91 [P3] CanSNPer2/CanSNPer2/modules/ModifyTree.py:405 — ValueError in update_annotations if row doesn't match expected format (sweep)"),
        ("P3 Medium", "- [ ] #92 [P3] CanSNPer2/CanSNPer2/modules/ModifyTree.py:458 — KeyError in update_genomes if genome or taxid is not in dictionaries (sweep)"),
        ("P3 Medium", "- [ ] #93 [P3] CanSNPer2/CanSNPer2/modules/ModifyTree.py:460 — KeyError in update_genomes if taxid is not in nodeDict (sweep)"),
        ("P3 Medium", "- [ ] #94 [P3] CanSNPer2/CanSNPer2/modules/ModifyTree.py:513 — ValueError in update_node_names if row doesn't match expected format (sweep)"),
        ("P3 Medium", "- [ ] #95 [P3] CanSNPer2/CanSNPer2/modules/NewickTree.py:230 — ZeroDivisionError in _confirm_path if dist is 0 (sweep)"),
        ("P3 Medium", "- [ ] #96 [P3] CanSNPer2/CanSNPer2/modules/ParseXMFA.py:68 — KeyError in reverse_complement if DNA has unknown characters (sweep)"),
        ("P3 Medium", "- [ ] #97 [P3] CanSNPer2/CanSNPer2/modules/ParseXMFA.py:102 — IndexError in find_snps if target sequence is shorter than reference (sweep)"),
        ("P3 Medium", "- [ ] #98 [P3] CanSNPer2/CanSNPer2/modules/ParseXMFA.py:136 — IndexError in read_sequence if seqLines is too short (sweep)"),
        ("P3 Medium", "- [ ] #99 [P3] CanSNPer2/CanSNPer2/modules/ParseXMFA.py:188 — IndexError in run if snp_positions is empty (sweep)"),
        ("P3 Medium", "- [ ] #100 [P3] CanSNPer2/CanSNPer2/modules/ReadCanSNP_tables.py:51 — IndexError in add_node_annotation if data is shorter than headers (sweep)"),
        ("P3 Medium", "- [ ] #101 [P3] CanSNPer2/CanSNPer2/modules/ReadCanSNP_tables.py:104 — IndexError in parse_snp_annotations if data has <2 elements (sweep)"),
        ("P3 Medium", "- [ ] #102 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomy.py:125 — IndexError in parse_taxonomy if nodes has only one element (sweep)"),
        ("P3 Medium", "- [ ] #103 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomy.py:84 — IndexError in read_nodes if data is shorter than 2 elements (sweep)"),
        ("P3 Medium", "- [ ] #104 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomy.py:95 — IndexError in read_nodes if data is shorter than 2 elements (sweep)"),
        ("P3 Medium", "- [ ] #105 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomyCanSNPer.py:137 — IndexError in read_nodes if data is empty (sweep)"),
        ("P3 Medium", "- [ ] #106 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomyCanSNPer.py:141 — IndexError in read_nodes if data has <2 elements (sweep)"),
        ("P3 Medium", "- [ ] #107 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomyCanSNPer.py:155 — KeyError in read_nodes if child/parent not in taxonomy (sweep)"),
        ("P3 Medium", "- [ ] #108 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomyCanSNPer.py:208 — IndexError in add_SNP if i is out of bounds (sweep)"),
        ("P3 Medium", "- [ ] #109 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomyCanSNPer.py:213 — KeyError in add_SNP if parent node not in taxonomy (sweep)"),
        ("P3 Medium", "- [ ] #110 [P3] CanSNPer2/CanSNPer2/modules/ReadTaxonomyCanSNPer.py:237 — KeyError in parse_taxonomy if parent name not in taxonomy (sweep)"),
        ("P3 Medium", "- [ ] #111 [P3] CanSNPer2/CanSNPer2/modules/WriteTaxonomy.py:25 — IndexError in __init__ if prefix doesn't have 2 elements (sweep)"),
        ("P3 Medium", "- [ ] #112 [P3] CanSNPer2/CanSNPer2/modules/WriteTaxonomy.py:117 — ValueError in unique_indexes if nodeDict is empty (sweep)"),
        ("P3 Medium", "- [ ] #113 [P3] CanSNPer2/CanSNPer2/modules/WriteTaxonomy.py:139 — KeyError in nodes if link[0] not in nodeDict (sweep)"),
        ("P3 Medium", "- [ ] #114 [P3] CanSNPer2/CanSNPer2/modules/WriteTaxonomy.py:150 — KeyError in nodes if link[1] not in nodeDict (sweep)"),
        ("P3 Medium", "- [ ] #115 [P3] CanSNPer2/CanSNPer2/modules/runMauve.py:46 — KeyError in RunMauve.__init__ if required kwargs are missing (sweep)"),
        ("P3 Medium", "- [ ] #116 [P3] CanSNPer2/CanSNPer2/modules/runMauve.py:81 — fragile command parsing with split(\" \") in run_mauve (sweep)"),
        ("P3 Medium", "- [ ] #117 [P3] CanSNPer2/CanSNPer2/modules/runMauve.py:142 — cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
        ("P3 Medium", "- [ ] #118 [P3] CanSNPer2/CanSNPer2/modules/runMauve.py:143 — cleanup may fail if tmpdir is not empty or contains directories (sweep)"),
    ]
    for priority, line in issues:
        append_issue(file_path, priority, line)
